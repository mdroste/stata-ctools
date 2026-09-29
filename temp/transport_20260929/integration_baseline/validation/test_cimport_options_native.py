"""ASan/UBSan coverage for delimiter normalization and unbounded CSV strings."""
import os
import platform
import subprocess
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char delimiters[40],collapse[8],literal[8],hex[16001],next[40],encoding[40];
static ST_retcode macro(char *key,char *out,ST_int cap) {
    const char *value="";
    if(!strcmp(key,"___cimport_delimiters"))value=delimiters;
    if(!strcmp(key,"___cimport_collapse"))value=collapse;
    if(!strcmp(key,"___cimport_asstring"))value=literal;
    if(!strcmp(key,"___cimport_encoding_name"))value=encoding;
    snprintf(out,(size_t)cap+1,"%s",value);return 0;
}
static ST_retcode save(char *key,char *value) {
    if(!strcmp(key,"_cimport_blobhex"))strcpy(hex,value);
    if(!strcmp(key,"_cimport_blobnext"))strcpy(next,value);
    return 0;
}
static ST_retcode display(char *s){(void)s;return 0;}
static ST_int poll(void){return 0;}
double ctools_timer_ms(void){return 0;}
#include "cimport/cimport_impl.c"
static void normalize(const char *input,const char *expected,int mode) {
    CImportContext ctx={0};ctx.file_data=(char*)input;ctx.file_size=strlen(input);ctx.bindquotes=mode;
    assert(!cimport_normalize_delimiters(&ctx));
    for(size_t i=0;i<strlen(expected);i++)assert(ctx.file_data[i]==(expected[i]=='@'?ctx.delimiter:expected[i]));
    assert(ctx.file_size==strlen(expected));free(ctx.converted_data);
}
int main(void) {
    mock.macuse=macro;mock.macresave=save;mock.spoutsml=display;mock.spouterr=display;mock.pollstd=poll;
    uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    strcpy(delimiters,"|;");strcpy(collapse,"1");
    normalize("||a;;\"b|c\";\n","@a@\"b|c\"@\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(collapse,"0");normalize("a|b;c\n","a@b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(delimiters,"||");strcpy(literal,"1");
    normalize("a|b||c\n","a|b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(delimiters,"\xce\xa9|");strcpy(literal,"0");
    normalize("a\xce\xa9" "b|c\n","a@b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    char *converted=NULL;size_t converted_size=0;
    const char utf32[]={0,0,0,'a',0,0,0,'b'};
    assert(!cimport_convert_named_to_utf8(utf32,sizeof(utf32),"UTF-32BE",&converted,&converted_size));
    assert(converted_size==2 && !strcmp(converted,"ab"));free(converted);
    assert(cimport_convert_named_to_utf8("x",1,"not-a-charset",&converted,&converted_size)==198);
    assert(!cimport_convert_named_to_utf8("\xff",1,"UTF-8",&converted,&converted_size));
    assert(!converted_size);free(converted);
    assert(!strcmp(cimport_detect_charset(NULL,0),"UTF-8"));
    assert(!strcmp(cimport_detect_charset("id,text\n1,first\n2,second\n",25),"ISO-8859-1"));
    const char utf32bom[]={ (char)0xff,(char)0xfe,0,0,'a',0,0,0 };
    assert(!strcmp(cimport_detect_charset(utf32bom,sizeof(utf32bom)),"UTF-32LE"));
    /* Prefixes and random byte buffers exercise every recognizer under ASan. */
    char sample[16000];uint32_t seed=17;
    for(size_t i=0;i<sizeof(sample);i++) {
        seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;sample[i]=(char)seed;
    }
    for(size_t n=0;n<=8000;n+=31)assert(*cimport_detect_charset(sample,n));
    for(size_t n=0;n<64;n++)assert(*cimport_detect_charset(sample,n));
    const char *before=cimport_detect_charset(sample,8000);
    memset(sample+8000,'x',8000);
    assert(!strcmp(before,cimport_detect_charset(sample,sizeof(sample))));
    double value;
    assert(!cimport_configure_locale("ar_EG","",""));
    assert(cimport_parse_unquoted_number("١٬٢٣٤٫٥٦",(int)strlen("١٬٢٣٤٫٥٦"),&value,SV_missval,'.',0));
    assert(fabs(value-1234.56)<1e-10);
    assert(!cimport_configure_locale("fa_IR","",""));
    assert(cimport_parse_unquoted_number("‎−۱۲۳۴٫۵۶",(int)strlen("‎−۱۲۳۴٫۵۶"),&value,SV_missval,'.',0));
    assert(fabs(value+1234.56)<1e-10);
    assert(cimport_configure_locale("bogus","","")==198);
    assert(!cimport_configure_locale("en_US","",""));
    assert(!cimport_parse_unquoted_number("+123",4,&value,SV_missval,'.',0));
    assert(cimport_parse_unquoted_number("∞",3,&value,SV_missval,'.',0) && value==SV_missval);
    assert(!cimport_configure_locale("","",""));
    /* Empty files must still validate charsets; conversion receives NULL/0. */
    delimiters[0]=collapse[0]=literal[0]=0;
    char emptyfile[]="/tmp/ctools-cimport-empty-XXXXXX";
    int fd=mkstemp(emptyfile);assert(fd>=0);close(fd);
    CImportContext *empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(empty && empty->num_columns==0 && empty->delimiter==',');
    assert(!strcmp(empty->reported_encoding,"UTF-8"));cimport_free_context(empty);
    strcpy(encoding,"latin1");
    empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(empty && !strcmp(empty->reported_encoding,"ISO-8859-1"));cimport_free_context(empty);
    strcpy(encoding,"bogus");
    empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(!empty && g_cimport_parse_status==198);encoding[0]=0;unlink(emptyfile);
    /* Ignore decoded NULs without truncating fields or mutating the mapping.
     * Exercise both the mapped ASCII path and the owned conversion buffer. */
    const char zeros[]="i\0d,name\n1\0,te\0st\n2,\"n\0ormal\"\n";
    char nullfile[]="/tmp/ctools-cimport-nuls-XXXXXX";
    fd=mkstemp(nullfile);assert(fd>=0);
    assert(write(fd,zeros,sizeof(zeros)-1)==sizeof(zeros)-1);close(fd);
    for(int explicit=0;explicit<2;explicit++) {
        strcpy(encoding,explicit?"UTF-8":"");
        CImportContext *clean=cimport_parse_csv(nullfile,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(clean && clean->num_columns==2);
        const char expected[]="id,name\n1,test\n2,\"normal\"\n";
        assert(clean->file_size==sizeof(expected)-1);
        assert(!memcmp(clean->file_data,expected,sizeof(expected)-1));
        cimport_free_context(clean);
    }
    encoding[0]=0;unlink(nullfile);
    CImportContext numeric={0};numeric.file_data="1\"\"2";numeric.strip_quotes_all=true;
    CImportFieldRef number={0,4|CIMPORT_FIELD_QUOTED_FLAG};
    assert(cimport_context_number(&numeric,&number,&value) && value==12);
    numeric.file_data="NA";number=(CImportFieldRef){0,2};
    assert(!cimport_context_number(&numeric,&number,&value));
    numeric.file_data="1\"\"2";number=(CImportFieldRef){0,4|CIMPORT_FIELD_QUOTED_FLAG};
    numeric.strip_quotes_all=false;numeric.retain_quotes=true;
    assert(!cimport_context_number(&numeric,&number,&value));
    CImportContext limit={0};limit.file_data="1,\"a\nb\nc\"\n";limit.file_size=strlen(limit.file_data);
    limit.bindquotes=CIMPORT_BINDQUOTES_STRICT;limit.max_quoted_rows=2;
    assert(cimport_check_quoted_rows(&limit)==5101);limit.max_quoted_rows=3;assert(!cimport_check_quoted_rows(&limit));
    CImportContext *ctx=calloc(1,sizeof(*ctx));pthread_mutex_init(&ctx->warning_mutex,NULL);
    ctx->num_columns=1;ctx->total_rows=1;ctx->num_chunks=1;ctx->num_threads=1;ctx->quote_char='"';
    ctx->columns=calloc(1,sizeof(*ctx->columns));ctx->columns[0].type=CIMPORT_COL_STRING;ctx->columns[0].max_strlen=20000;
    ctx->columns[0].use_strl=true;
    ctx->chunks=calloc(1,sizeof(*ctx->chunks));ctx->chunks[0].num_rows=1;
    ctx->chunks[0].rows=calloc(1,sizeof(CImportParsedRow*));
    CImportParsedRow *row=ctools_arena_alloc(&ctx->chunks[0].arena,cimport_row_size(1));
    CImportFieldRef long_field={0,20000};
    assert(row && cimport_fill_row(row,&long_field,1));cimport_row_values(row)[0]=CIMPORT_NO_VALUE;
    ctx->chunks[0].rows[0]=row;
    ctx->file_data=ctx->converted_data=malloc(20001);ctx->file_size=20000;memset(ctx->file_data,'x',20000);ctx->file_data[20000]=0;
    atomic_init(&ctx->error_code,0);g_cimport_ctx=ctx;
    assert(!cimport_blob(1,1,0));assert(strlen(hex)==16000 && !strcmp(next,"8000"));
    assert(!cimport_blob(1,1,16000));assert(strlen(hex)==8000 && !strcmp(next,"0"));
    assert(strlen(ctx->col_cache[0].string_data[0])==20000);
    assert(cimport_blob(1,1,20001)==198);assert(cimport_blob(2,1,0)==198);
    cimport_free_context(ctx);g_cimport_ctx=NULL;
    ctools_destroy_global_pool();return 0;
}
'''

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-cimport-options-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'; source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined',
                 '-fno-omit-frame-pointer', '-pthread', '-I', str(ROOT / 'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                      '-Wl,-undefined,dynamic_lookup', '-liconv']
        else:
            flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                      '-fdata-sections', '-Wl,--gc-sections', '-lm']
        modules=['ctools_types.c','ctools_arena.c','ctools_threads.c',
                 'cimport/cimport_parse.c','cimport/cimport_mmap.c','cimport/cimport_encoding.c']
        subprocess.run([os.environ.get('CC', 'cc'), *flags, str(source),
                        *[str(ROOT/'src'/m) for m in modules], '-o', str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60,env={**os.environ,'UBSAN_OPTIONS':'halt_on_error=1'})
        print('PASS delimiter sets/literals/Unicode, quote limits, long CSV cache/blob bounds (ASan/UBSan)')

if __name__ == '__main__':
    main()

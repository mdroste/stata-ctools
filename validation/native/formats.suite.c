/* Formats tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_io_formats_native_source)
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
#include "io/cio_context.c"
#include "io/cio_readstat.c"
#include "io/cio_dbf.c"
#include "io/cio_metadata.c"
#include "io/cio_excel.c"
#include "io/cio_xls.c"
#include "io/cio_xlsx.c"
#include "io/cio_xport.c"
#include "io/cio_xpf.c"
#include "io/cio_shp.c"
#include "formats/excel_address.c"
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
    assert(cio_scan(&cache, kind)==0);assert(cache.nvar==columns && cache.nobs==rows);
    if(!strcmp(kind,"spss") || !strcmp(kind,"sas")){
        assert(strlen(cache.cols[2].strings[0])==2100);
        assert(cache.cols[2].width==2100);
        char type[40];assert(!strcmp(cio_storage(&cache, &cache.cols[2],type,sizeof(type)),"strL"));
        assert(!strcmp(cache.cols[2].format,"%9s"));
        assert(cio_blob_info(&cache, 3,1,0)==0);assert(cio_blob_info(&cache, 3,1,2101)==198);
        assert(cio_load(&cache)==0);store_error=1;assert(cio_load(&cache)==459);store_error=0;
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
    assert(!cio_scan(&cache, "sasxport8"));assert(cache.nvar==2 && cache.nobs==3);
    char type[40];assert(!strcmp(cio_storage(&cache, &cache.cols[0],type,sizeof(type)),"byte"));
    assert(cache.cols[1].width==3);cio_cleanup();raw_storage=1;
    assert(!cio_scan(&cache, "sasxport8"));assert(cache.cols[1].width==12);
    assert(!strcmp(cio_storage(&cache, &cache.cols[0],type,sizeof(type)),"double"));
    assert(!strcmp(cache.cols[1].format,"%12s"));cio_cleanup();raw_storage=0;
    for(int i=0;i<2;i++) {
        snprintf(filename,sizeof(filename),"%s/display_formats.%s",argv[1],i ? "xlsx" : "xls");
        assert(!cio_scan(&cache, i ? "xlsx" : "xls"));assert(cache.nvar==36 && cache.nobs==4);
        assert(!strcmp(cache.cols[27].format,"%10.0g"));
        assert(cache.cols[27].numbers[0]==44300.5123456);
        assert(cache.cols[35].numbers[0]==44300.5123456);
        assert(!strcmp(cache.cols[4].format,"%tchh:MM_AM"));
        cio_cleanup();excel_strings=1;
        assert(!cio_scan(&cache, i ? "xlsx" : "xls"));
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
        assert(cio_scan(&cache, kinds[j])!=0);cio_cleanup();
    }
    /* Truncated SAS pages and ZSAV blocks must not leave partially allocated
       columns or dangling codec state across repeated import attempts. */
    for(int file=0;file<2;file++){
        char original[4096];snprintf(original,sizeof(original),"%s/%s",argv[1],file?"fixture.zsav":"long.sas7bdat");
        FILE *input=fopen(original,"rb");assert(input);unsigned char bytes[16384];size_t n=fread(bytes,1,sizeof(bytes),input);fclose(input);
        for(size_t length=0;length<n;length+=97){
            snprintf(filename,sizeof(filename),"%s/truncated.bin",argv[1]);FILE *out=fopen(filename,"wb");assert(out);assert(fwrite(bytes,1,length,out)==length);assert(!fclose(out));
            (void)cio_scan(&cache, file?"spss":"sas");cio_cleanup();
        }
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_io_parity_native_source)
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_boolean missing(double x) { return x >= mock.missval; }
#include "cexport/cexport_format.c"
#include "cimport/cimport_parse.c"
int main(void) {
    uint64_t bits = UINT64_C(0x7fe0000000000000);
    memcpy(&mock.missval, &bits, 8);
    mock.ismissing = missing;
    char buffer[8192], source[2046];
    for (int i=0; i<2045; i++) source[i] = i % 2 ? 'x' : '"';
    source[2045] = 0;
    int n = cexport_write_quoted_string(source, buffer, sizeof(buffer));
    assert(n == 2045 + 1023 + 2);
    assert(buffer[0] == '"' && buffer[n-1] == '"');
    assert(cexport_write_quoted_string(source, buffer, 256) == -1);
    for (int i=1; i<=26; i++) {
        uint64_t encoded = bits + ((uint64_t)i << 40);
        double x, parsed;
        memcpy(&x, &encoded, 8);
        assert(cexport_double_to_str(x, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) == 2);
        assert(buffer[0] == '.' && buffer[1] == 'a'+i-1);
        assert(cimport_field_looks_numeric_sep(buffer, 2, '.', 0));
        assert(cimport_parse_number(buffer, 2, &parsed, mock.missval, '.', 0));
        assert(parsed == x);
    }
    assert(!cimport_field_looks_numeric_sep("NA", 2, '.', 0));
    assert(!cimport_field_looks_numeric_sep("NaN", 3, '.', 0));
    double number;
    assert(cimport_parse_number("\"1.25\"", 6, &number, mock.missval, '.', 0) && number == 1.25);
    const char *raw = "1e-2";
    CImportFieldRef field = {0, 4};
    bool integer;
    assert(cimport_analyze_numeric_with_sep(raw, &field, '"', '.', 0, &number, &integer));
    assert(!integer && number == 0.01);
    assert(cexport_double_to_str(1e100, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer,"1.0000000000e+100"));
    assert(cexport_double_to_str(1e-5, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer,".00001"));
    assert(cexport_double_to_str(-27915931.612694785, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer, "-27915931.61269479"));
    assert(cexport_double_to_str(136948012240.34935, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer, "136948012240.3494"));
    /* A long quoted field must also survive the mixed-type row writer. */
    stata_variable var = {0}; char *strings[] = {source};
    var.type = STATA_TYPE_STRING; var.data.str = strings;
    vartype_t type = VARTYPE_STRING;
    cexport_context ctx = {0};
    ctx.filtered.data.nvars = 1; ctx.filtered.data.vars = &var;
    ctx.vartypes = &type; ctx.delimiter = ','; ctx.quote_if_needed = true;
    ctx.line_ending[0] = '\n'; ctx.line_ending_len = 1;
    assert(cexport_format_row(&ctx, 0, buffer, sizeof(buffer)) == n + 1);
    assert(buffer[n] == '\n');
    assert(cexport_format_row(&ctx, 0, buffer, 256) == -1);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_xls_edit_native_source)
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
#endif

#if defined(CTOOLS_CASE_xlsx_edit_native_source)
#include <assert.h>
#include <string.h>
#include "stplugin.h"
ST_plugin *_stata_;
#include "cexport/cexport_xlsx.c"
#include "ctools_text.c"
#include "formats/excel_address.c"
int main(void) {
    const char *old = "<x:worksheet xmlns:x=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
        "<!-- > <row> --> <x:dimension ref='A1:F9'/><x:sheetData>"
        "<x:row r='2' ht='30' customHeight='1'><x:c r='A2' s='8'><x:f>SUM(A1)&gt;0</x:f><x:v>9</x:v></x:c>"
        "<x:c r='B2' s='7'><x:v>3</x:v></x:c><x:c r='F2'><x:is><x:t>old</x:t></x:is></x:c></x:row>"
        "<x:row r='9'><x:c r='F9'><x:v>12</x:v></x:c></x:row></x:sheetData>"
        "<x:mergeCells count='1'><x:mergeCell ref='G10:H10'/></x:mergeCells></x:worksheet>";
    const char *fresh="<worksheet><sheetData><row r='2'><c r='B2' s='1'><v>20</v></c><c r='C2'/></row>"
        "<row r='3'><c r='B3'><is><t>new &amp; text</t></is></c></row></sheetData></worksheet>";
    char *got=xedit_worksheet(old,fresh,1,20);
    assert(got && strstr(got,"s=\"7\"") && strstr(got,"s=\"20\""));
    assert(strstr(got,"<x:c r='A2' s='8'><x:f>SUM(A1)&gt;0</x:f><x:v>9</x:v></x:c>"));
    assert(strstr(got,"<x:c r='F9'><x:v>12</x:v></x:c>"));
    assert(strstr(got,"ht='30' customHeight='1'"));
    assert(strstr(got,"ref=\"A2:F9\"") && strstr(got,"ref='G10:H10'"));
    assert(!strstr(got,"<x:v>3</x:v>")); free(got);
    got=xedit_worksheet(old,fresh,0,20);assert(got && strstr(got,"s=\"21\""));free(got);
    got=xedit_worksheet("<worksheet><sheetData/></worksheet>",fresh,0,0);
    assert(got && strstr(got,"new &amp; text"));free(got);
    xedit_node n;
    const char *entities="<sheet name='A &amp; B &#x3a9;'/>";
    assert(xedit_find(entities,entities+strlen(entities),"sheet",&n));
    got=xedit_attr(&n,"name");assert(got && !strcmp(got,"A & B \xce\xa9"));free(got);
    const char *styles="<styleSheet><numFmts count='1'><numFmt numFmtId='181' formatCode='custom'/></numFmts>"
        "<cellXfs count='1'><xf numFmtId='181' fontId='2' fillId='3'/></cellXfs></styleSheet>";
    unsigned base=99;got=xedit_styles((char*)styles,&base);
    assert(got && base==1 && strstr(got,"numFmtId=\"182\"") && strstr(got,"count=\"4\""));
    assert(strstr(got,"<xf numFmtId='181' fontId='2' fillId='3'/>"));free(got);
    styles="<styleSheet><cellXfs count='0'/></styleSheet>";
    got=xedit_styles((char*)styles,&base);assert(got && base==0 && strstr(got,"count=\"3\""));free(got);
    got=xedit_translate_formula("SUM(A2,$A2,A$2,$A$2,'A2 sheet'!B3,LOG10(A2),\"A2\",Table1[A2],A2:B4)",1,2);
    assert(got && !strcmp(got,"SUM(B4,$A4,B$2,$A$2,'A2 sheet'!C5,LOG10(B4),\"A2\",Table1[A2],B4:C6)"));free(got);
    got=xedit_translate_formula("A1+$A1+A$1+$A$1",-1,-1);
    assert(got && !strcmp(got,"#REF!+#REF!+#REF!+$A$1"));free(got);
    got=xedit_translate_formula("SUM($1:2,A:$B,'A1:A3'!C2,A1:A3!B2)",1,2);
    assert(got && !strcmp(got,"SUM($1:2,A:$B,'A1:A3'!D4,A1:A3!C4)"));free(got);
    const char *shared="<worksheet><sheetData><row r='2'><c r='C2'><f t='shared' ref='C2:C4' si='0'>A2+B2</f><v>5</v></c></row>"
        "<row r='3'><c r='C3'><f t='shared' si='0'/><v>7</v></c></row><row r='4'><c r='C4'><f t='shared' si='0'/><v>9</v></c></row></sheetData></worksheet>";
    const char *overwritten="<worksheet><sheetData><row r='2'><c r='C2'><v>99</v></c></row>"
        "<row r='3'><c r='C3'><f t='shared' si='0'/><v>7</v></c></row><row r='4'><c r='C4'><f t='shared' si='0'/><v>9</v></c></row></sheetData></worksheet>";
    got=xedit_shared_anchors(shared,overwritten);
    assert(got && strstr(got,"<f t='shared' ref='C2:C4' si='0'>A3+B3</f>"));
    assert(strstr(got,"<v>99</v>") && strstr(got,"<c r='C4'><f t='shared' si='0'/>"));free(got);
    char chain[100];
    got=xedit_remove_chain("<Relationships><Relationship Id='one' Type='url/worksheet' Target='s.xml'/><Relationship Id='two' Type='url/calcChain' Target='calcChain.xml'/></Relationships>","Relationship","Type","/calcChain",chain,sizeof(chain));
    assert(got && !strcmp(chain,"xl/calcChain.xml") && strstr(got,"Id='one'") && !strstr(got,"Id='two'"));free(got);
    got=xedit_recalculate("<x:workbook xmlns:x='urn:x'><x:calcPr fullCalcOnLoad='0' calcId='42'/></x:workbook>");
    assert(got && strstr(got,"fullCalcOnLoad=\"1\"") && strstr(got,"calcId='42'"));free(got);
    got=xedit_recalculate("<workbook/>");
    assert(got && strstr(got,"fullCalcOnLoad=\"1\"") && strstr(got,"</workbook>"));free(got);
    const char *lex="SUM($1:2,A:$B,'A1:A3'!C2,A1:A3!B2,Table1[[A2]:[B3]],\"A2\"\"B3\",a2)";
    for(size_t i=0;i<=strlen(lex);i++) {
        char *truncated=malloc(i+1);memcpy(truncated,lex,i);truncated[i]=0;
        got=xedit_translate_formula(truncated,1,1);assert(got);free(got);free(truncated);
    }
    /* Truncated tags must be bounded: no reads past supplied lengths. */
    for(size_t i=0;i<strlen(old);i++) {
        char *truncated=malloc(i+1);memcpy(truncated,old,i);truncated[i]=0;
        xedit_find(truncated,truncated+i,"sheetData",&n);
        got=xedit_worksheet(truncated,fresh,1,0);free(got);free(truncated);
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_xlsx_strl_native_source)
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
#include "ctools_text.c"
#include "formats/excel_address.c"
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
    assert(!ctools_load_text_data(&fd,3,(ctools_text_policy){0, 0}));
    assert(strlen(fd.data.vars[0].data.str[0])==100000);
    ctools_filtered_data_free(&fd);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_oct21_xlsx)
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
double ctools_timer_seconds(void) { return 0; }
#include "formats/xlsx_zip.c"
#include "cimport/cimport_xlsx.c"

static void scan(const char *xml, int expected)
{
    XLSXContext ctx;
    xlsx_context_init(&ctx);
    /* A range must not hide corrupt XML outside the selection. */
    ctx.cell_range.end_row = 1;
    assert(xlsx_scan_buffer(&ctx, xml, strlen(xml)) == expected);
    xlsx_context_free(&ctx);
}
int main(int argc, char **argv)
{
    assert(argc == 2);
    uint64_t missing = UINT64_C(0x7fe0000000000000);
    memcpy(&mock.missval, &missing, sizeof(missing));
    ctools_set_max_threads(4);
    const char *files[] = {"valid.xlsx", "stored.xlsx", "bitflip.xlsx", "length.xlsx", "crc.xlsx"};
    for (int i = 0; i < 5; i++) {
        char path[4096]; snprintf(path, sizeof(path), "%s/%s", argv[1], files[i]);
        xlsx_zip_archive *zip = xlsx_zip_open_file(path); assert(zip);
        size_t index = xlsx_zip_locate_file(zip, "sheet.xml");
        for (int direct = 0; direct < 2; direct++) {
            size_t size = 0;
            void *bytes = direct ? xlsx_zip_extract_direct(zip, index, &size)
                                 : xlsx_zip_extract_fast(zip, index, &size);
            assert((bytes != NULL) == (i < 2));
            if (bytes) { assert(size == 10); assert(!memcmp(bytes, "<v>123</v>", 10)); }
            free(bytes);
        }
        xlsx_zip_close(zip);
    }
    const char *good = "<worksheet><sheetData><row r=\"1\"><c r=\"A1\"><v>123</v></c></row></sheetData></worksheet>";
    scan(good, 0);
    /* Every premature EOF, including between complete rows, must fail. */
    char truncated[256]; strcpy(truncated, good);
    for (size_t i = 0; i < strlen(good); i++) {
        memcpy(truncated, good, i); truncated[i] = 0; scan(truncated, 610);
    }
    scan("<worksheet><sheetData/></worksheet>", 0);
    scan("<?xml version='1.0'?><s:worksheet xmlns:s='n'><!-- ok --><s:sheetData><s:row r='1'><s:c/></s:row></s:sheetData></s:worksheet>", 0);
    scan("<worksheet><sheetData><row><c><v>1</c></v></row></sheetData></worksheet>", 610);
    scan("<worksheet><sheetData/><extra x='unterminated></worksheet>", 610);
    scan("<worksheet><sheetData/></worksheet><!-- unfinished", 610);
    /* Real parallel row partitions and an early cellrange termination. */
    size_t capacity = 500000, used = 0;
    char *large = malloc(capacity); assert(large);
    used += (size_t)snprintf(large + used, capacity - used, "<worksheet><sheetData>");
    for (int row = 1; row <= 4000; row++)
        used += (size_t)snprintf(large + used, capacity - used,
            "<row r=\"%d\"><c r=\"A%d\"><v>%d</v></c></row>", row, row, row);
    used += (size_t)snprintf(large + used, capacity - used, "</sheetData></worksheet>");
    assert(used < capacity && used > XLSX_PARALLEL_MIN_BYTES);
    for (int range = 0; range < 2; range++) {
        XLSXContext ctx; xlsx_context_init(&ctx);
        if (range) ctx.cell_range.end_row = 1;
        assert(!xlsx_scan_buffer(&ctx, large, used));
        assert(ctx.cm_numeric[0][0] == 1);
        if (!range) assert(ctx.cm_numeric[0][3999] == 4000);
        xlsx_context_free(&ctx);
    }
    large[used - 1] = 0; scan(large, 610); free(large);
    ctools_destroy_global_pool();
    puts("PASS ZIP integrity and complete worksheet XML");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_oct21_output)
#include "ctools_file.h"
#include <assert.h>
#include <unistd.h>
static const char *fail_destination;
static int force_copy;
#ifdef __APPLE__
static int test_rename_excl(const char *from, const char *to, unsigned flags)
{
    if (fail_destination && !strcmp(to, fail_destination)) {
        fail_destination = NULL; errno = EIO; return -1;
    }
    if (force_copy) { errno = ENOTSUP; return -1; }
    return renamex_np(from, to, flags);
}
#define renamex_np test_rename_excl
#elif defined(__linux__)
#include <sys/syscall.h>
#include <stdarg.h>
/* Match syscall's variadic declaration when the implementation is included. */
static long test_syscall(long number, ...)
{
    assert(number == SYS_renameat2);
    va_list args;
    va_start(args, number);
    int dir1 = va_arg(args, int);
    const char *from = va_arg(args, const char *);
    int dir2 = va_arg(args, int);
    const char *to = va_arg(args, const char *);
    int flags = va_arg(args, int);
    va_end(args);
    if (fail_destination && !strcmp(to, fail_destination)) {
        fail_destination = NULL; errno = EIO; return -1;
    }
    if (force_copy) { errno = ENOTSUP; return -1; }
    return syscall(number, dir1, from, dir2, to, flags);
}
#define syscall test_syscall
#endif
static int test_link(const char *from, const char *to)
{
    if (force_copy) { errno = ENOTSUP; return -1; }
    return link(from, to);
}
#define link test_link
#include "formats/ctools_output.c"

static void write_file(const char *path, const char *text)
{
    FILE *fp = fopen(path, "wb"); assert(fp);
    assert(fwrite(text, 1, strlen(text), fp) == strlen(text)); assert(!fclose(fp));
}
static void contents(const char *path, const char *text)
{
    char buf[100] = {0}; FILE *fp = fopen(path, "rb"); assert(fp);
    fread(buf, 1, sizeof(buf)-1, fp); fclose(fp); assert(!strcmp(buf, text));
}
static void prepare(ctools_output *out, const char *path, int replace)
{
    assert(!ctools_output_prepare(out, path, replace)); write_file(out->temporary, "new");
}
int main(int argc, char **argv)
{
    assert(argc == 2); assert(!chdir(argv[1])); umask(0027);
    ctools_output a, b; ctools_output *pair[] = {&a, &b};
    prepare(&a, "mode", 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    struct stat st; assert(!stat("mode", &st)); assert((st.st_mode & 0777) == 0640);
    assert(!chmod("mode", 0664)); prepare(&a, "mode", 1);
    assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    assert(!stat("mode", &st)); assert((st.st_mode & 0777) == 0664);
    char name[251]; memset(name, 'x', 250); name[250] = 0;
    prepare(&a, name, 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a); contents(name, "new");
    force_copy = 1; prepare(&a, "fat", 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    contents("fat", "new"); prepare(&a, "fat", 0); assert(ctools_output_commit(&a) == 602); ctools_output_cleanup(&a);
    force_copy = 0;
    write_file("b", "old-b"); prepare(&a, "a", 0); prepare(&b, "b", 0);
    assert(ctools_output_commit_group(pair, 2) == 602); assert(access("a", F_OK)); contents("b", "old-b");
    ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    for (int existing = 0; existing < 2; existing++) {
        if (existing) write_file("a", "old-a");
        prepare(&a, "a", existing); prepare(&b, "b", 1);
        fail_destination = "b";
        assert(ctools_output_commit_group(pair, 2) == 603);
        if (existing) contents("a", "old-a"); else assert(access("a", F_OK));
        contents("b", "old-b"); ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    }
    prepare(&a, "a", 1); prepare(&b, "b", 1); assert(!ctools_output_commit_group(pair, 2));
    contents("a", "new"); contents("b", "new"); ctools_output_cleanup(&a); ctools_output_cleanup(&b);
#ifdef __APPLE__
    write_file("acl", "old");
    assert(!system("chmod +a 'everyone allow read' acl"));
    acl_t acl = acl_get_file("acl", ACL_TYPE_EXTENDED); assert(acl);
    char *before = acl_to_text(acl, NULL); assert(before); acl_free(acl);
    prepare(&a, "acl", 1); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    acl = acl_get_file("acl", ACL_TYPE_EXTENDED); assert(acl);
    char *after = acl_to_text(acl, NULL); assert(after); assert(!strcmp(before, after));
    acl_free(before); acl_free(after); acl_free(acl);
#endif
    write_file("alias", "original"); prepare(&a, "alias", 1); prepare(&b, "./alias", 1);
    assert(ctools_output_commit_group(pair, 2) == 198); contents("alias", "original");
    ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    puts("PASS output rollback, permissions/ACLs, long names, and no-hardlink fallback");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_oct21_windows)
/* Run with the release compiler on a real Windows worker, without Stata. */
#include "ctools_file.h"
#include "formats/ctools_output.h"
#include <assert.h>
#include <locale.h>

static void roundtrip(const char *path)
{
    ctools_output out;
    assert(!ctools_output_prepare(&out, path, false));
    FILE *fp = ctools_fopen(out.temporary, "wb"); assert(fp);
    assert(fwrite("123", 1, 3, fp) == 3); assert(!fclose(fp));
    assert(!ctools_output_commit(&out)); ctools_output_cleanup(&out);
    fp = ctools_fopen(path, "rb"); assert(fp); char buf[4] = {0};
    assert(fread(buf, 1, 3, fp) == 3); assert(!strcmp(buf, "123")); fclose(fp);
    assert(!ctools_output_prepare(&out, path, false));
    assert(ctools_output_commit(&out) == 602); ctools_output_cleanup(&out);
    assert(!ctools_output_prepare(&out, path, true));
    fp = ctools_fopen(out.temporary, "wb"); assert(fp); assert(!fclose(fp));
    assert(!ctools_output_commit(&out)); ctools_output_cleanup(&out);
    HANDLE h = ctools_create_file(path, GENERIC_READ, FILE_SHARE_READ, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL);
    assert(h != INVALID_HANDLE_VALUE); CloseHandle(h);
    assert(!ctools_remove(path));
}
int main(int argc, char **argv)
{
    assert(argc == 2);
    /* Force a legacy CRT locale, independently of the machine's UI language. */
    assert(setlocale(LC_ALL, "C"));
    char path[4096], parent[4096];
    snprintf(parent, sizeof(parent), "%s/\xc3\xa9-\xe6\xb8\xac-\xf0\x9f\x98\x80", argv[1]);
    wchar_t *wide = ctools_wpath(parent); assert(wide);
    assert(CreateDirectoryW(wide, NULL)); free(wide);
    snprintf(path, sizeof(path), "%s/file.csv", parent); roundtrip(path);
    char extended[4096];
    wide = ctools_wpath(path); assert(wide);
    assert(WideCharToMultiByte(CP_UTF8, 0, wide, -1, extended, sizeof(extended), NULL, NULL));
    free(wide); roundtrip(extended);
    /* Nested paths longer than MAX_PATH, with ordinary UTF-8 path syntax. */
    size_t lengths[8];
    for (int i = 0; i < 8; i++) {
        lengths[i] = strlen(parent);
        strcat(parent, "/abcdefghijklmnopqrstuvwxyz0123456789");
        wide = ctools_wpath(parent); assert(wide); assert(CreateDirectoryW(wide, NULL)); free(wide);
    }
    snprintf(path, sizeof(path), "%s/file.xlsx", parent); roundtrip(path);
    FILE *fp = ctools_fopen(path, "w+b"); assert(fp);
    int64_t offset = INT64_C(2147483648) + 127;
    assert(!ctools_fseek(fp, offset, SEEK_SET)); assert(ctools_ftell(fp) == offset);
    fclose(fp); assert(!ctools_remove(path));
    for (int i = 7; i >= 0; i--) {
        wide = ctools_wpath(parent); assert(RemoveDirectoryW(wide)); free(wide);
        parent[lengths[i]] = 0;
    }
    wide = ctools_wpath(parent); assert(RemoveDirectoryW(wide)); free(wide);
    assert(!ctools_wpath("invalid-\xff"));
    puts("PASS Windows Unicode/long paths, replacement, cleanup, and 64-bit offsets");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_oct22_xml)
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
double ctools_timer_seconds(void) { return 0; }
#include "formats/xlsx_zip.c"
#include "cimport/cimport_xlsx.c"
static void text_case(const char *xml, const char *expected) {
    XLSXContext ctx;
    xlsx_context_init(&ctx);
    assert(xlsx_scan_buffer(&ctx,xml,strlen(xml)) == 0);
    assert(ctx.cm_types[0][0] == XLSX_CELL_STRING);
    uintptr_t p; memcpy(&p,&ctx.cm_numeric[0][0],sizeof(p));
    assert(!strcmp((const char *)p,expected));
    xlsx_context_free(&ctx);
}
int main(void) {
    uint64_t missing=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&missing,8);
    const char *sst="<sst><si><t>first</t></si><si><t>second";
    XLSXContext ctx;xlsx_context_init(&ctx);
    assert(xlsx_parse_shared_strings_buf(&ctx,sst,strlen(sst)) != 0);
    assert(ctx.num_shared_strings == 0);
    xlsx_context_free(&ctx);
    const char *styles="<styleSheet><cellXfs><xf numFmtId='14'/>";
    xlsx_context_init(&ctx);
    assert(xlsx_parse_styles_buf(&ctx,styles,strlen(styles)) != 0);
    xlsx_context_free(&ctx);
    text_case("<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t><![CDATA[customer &amp; name]]></t></is></c></row></sheetData></worksheet>","customer &amp; name");
    text_case("<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t>before<!-- > < -->after</t></is></c></row></sheetData></worksheet>","beforeafter");
    text_case("<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t>a_x000D_b_x005F_x0041_</t></is></c></row></sheetData></worksheet>","a_x000D_b_x005F_x0041_");
    const char *sst2="<sst><si><t><![CDATA[customer &amp; name]]></t></si></sst>";
    xlsx_context_init(&ctx);assert(xlsx_parse_shared_strings_buf(&ctx,sst2,strlen(sst2))==0);
    assert(!strcmp(ctx.shared_strings[0],"customer &amp; name"));xlsx_context_free(&ctx);
    char escaped[]="_xD83D__xDE00_";
    assert(xlsx_decode_ooxml(escaped,strlen(escaped))==4);
    assert(!strcmp(escaped,"\xf0\x9f\x98\x80"));
    ctools_destroy_global_pool();puts("oct22 XML passed");return 0;
}
#endif

#if defined(CTOOLS_CASE_oct22_export)
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
#include "cexport/cexport_xlsx.c"
#include "ctools_text.c"
#include "formats/excel_address.c"
int main(void)
{
    char escaped[128];
    assert(xlsx_escape_xml("a\001b\r_x0041_",escaped,sizeof(escaped))>0);
    assert(!strcmp(escaped,"a_x0001_b_x000D__x005F_x0041_"));
    XLSXExportContext ctx;
    xlsx_context_init(&ctx);
    assert(xlsx_collect_string(&ctx,"ordinary & text")==0);
    assert(ctx.num_shared_strings==0);
    assert(xlsx_collect_string(&ctx,"a\001b")==0);
    assert(xlsx_collect_string(&ctx,"a\001b")==0);
    assert(ctx.num_shared_strings==1);
    StringBuffer cell;
    assert(strbuf_init(&cell,128));
    ctx.shared_offset=7;
    assert(xlsx_shared_cell(&cell,&ctx,"<c r=\"A",1,"a\001b"));
    assert(!strcmp(cell.data,"<c r=\"A1\" t=\"s\"><v>7</v></c>\n"));
    free(cell.data);
    char *long_text=malloc(32768);
    assert(long_text);
    memset(long_text,1,32767);long_text[32767]=0;
    assert(xlsx_collect_string(&ctx,long_text)==0);
    free(long_text);
    char *items=xlsx_shared_items(&ctx);
    assert(items);
    assert(strlen(items)>32767*7);
    assert(!strcmp(items+strlen(items)-9,"</t></si>"));
    free(items);
    xlsx_context_free(&ctx);
    puts("oct22 OOXML export and long shared strings passed");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_styles)
#include "stplugin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
ST_plugin mock;
ST_plugin *_stata_ = &mock;
#include "alloc.h"
#include "cimport/cimport_xlsx.c"
#include "formats/xlsx_xml.c"
#include "formats/excel_address.c"
int main(int argc, char **argv) {
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    XLSXContext ctx = {0};
    const char *s = "<styleSheet><cellXfs count=\"2\"><xf numFmtId=\"0\"/><xf "
                    "numFmtId=\"14\"/></cellXfs></styleSheet>";
    int rc = xlsx_parse_styles_buf(&ctx, s, strlen(s));
    assert(fail_at ? rc == 920 : (rc == 0 && ctx.num_styles >= 2 && ctx.date_styles[1] == 1));
    printf("failure=%d rc=%d allocations=%d styles=%d date=%d\n", fail_at, rc, calls,
           ctx.num_styles, ctx.num_styles ? ctx.date_styles[1] : -1);
    free(ctx.date_styles);
    free(ctx.excel_formats);
    assert(alive == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_io_perf_crc)
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
#endif

#if defined(CTOOLS_CASE_io_perf_xml)
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
        if (type==XLSX_CELL_SHARED_STRING) {
            const char *sst="<sst><si><t>zero</t></si><si><t>one</t></si><si><t>two</t></si><si><t>three</t></si></sst>";
            assert(xlsx_parse_shared_strings_buf(&x,sst,strlen(sst))==0);
        }
        assert(xlsx_cm_reserve(&x,1,1));
        xlsx_scan_local loc; xlsx_scan_local_init(&loc);
        loc.arena=&x.parse_arena;loc.grow=&x;
        const char *end=input+length;
        const char *gt=length ? xlsx_markup_end(input,end) : NULL;
        int prev=-1;
        const char *after=gt ? xlsx_scan_cell(&x,&loc,input+2,gt,end,1,&prev,false) : NULL;
        if (length==n) {
            assert(after==end && !loc.failed && !loc.invalid);
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
#endif

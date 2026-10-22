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

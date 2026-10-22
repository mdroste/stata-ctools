#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
double ctools_timer_seconds(void) { return 0; }
#include "formats/xlsx_zip.c"
#include "cimport/cimport_xlsx.c"

static void scan_text(const char *tag, const char *xml) {
    XLSXContext ctx;
    xlsx_context_init(&ctx);
    int rc=xlsx_scan_buffer(&ctx,xml,strlen(xml));
    const char *text="";
    if(ctx.cm_active && ctx.cm_num_rows && ctx.cm_num_cols && ctx.cm_types[0][0]==XLSX_CELL_STRING) {
        uintptr_t p;memcpy(&p,&ctx.cm_numeric[0][0],sizeof(p));text=(const char*)p;
    }
    printf("%s: rc=%d text='%s' cells=%zu\n",tag,rc,text,ctx.cells_stored);
    xlsx_context_free(&ctx);
}
int main(void) {
    uint64_t missing=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&missing,8);
    const char *sst="<sst><si><t>first</t></si><si><t>second";
    XLSXContext ctx;xlsx_context_init(&ctx);
    int rc=xlsx_parse_shared_strings_buf(&ctx,sst,strlen(sst));
    printf("Truncated shared strings: rc=%d strings=%u first='%s'\n",rc,ctx.num_shared_strings,ctx.shared_strings[0]);
    xlsx_context_free(&ctx);
    scan_text("CDATA","<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t><![CDATA[customer & name]]></t></is></c></row></sheetData></worksheet>");
    scan_text("Comment text","<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t>before<!-- comment -->after</t></is></c></row></sheetData></worksheet>");
    scan_text("OOXML escaped text","<worksheet><sheetData><row r='1'><c r='A1' t='inlineStr'><is><t>a_x000D_b</t></is></c></row></sheetData></worksheet>");
    const char *sst2="<sst><si><t><![CDATA[customer & name]]></t></si></sst>";
    xlsx_context_init(&ctx);rc=xlsx_parse_shared_strings_buf(&ctx,sst2,strlen(sst2));
    printf("Shared CDATA: rc=%d text='%s'\n",rc,ctx.shared_strings[0]);xlsx_context_free(&ctx);
    ctools_destroy_global_pool();return 0;
}

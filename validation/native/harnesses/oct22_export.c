#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
#include "cexport/cexport_xlsx.c"
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

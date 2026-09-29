#include "stplugin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
ST_plugin mock;
ST_plugin *_stata_ = &mock;
#include "alloc.h"
#include "cimport/cimport_xlsx.c"
#include "cimport/cimport_xlsx_xml.c"
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

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

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

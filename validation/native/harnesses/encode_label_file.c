#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ctools_hash.c"

int main(void)
{
    char path[] = "/tmp/ctools-labels-XXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);
    close(fd);
    char longest[2046];
    memset(longest, 'x', 2045);
    longest[2045] = '\0';
    const char *texts[] = {"", "\r\n\t", "`macro' $global \"quoted\"", "é東京", longest};
    const int codes[] = {-10, 0, 1, 7, 65536};
    assert(ctools_label_write_file(texts, codes, 5, path) == 0);
    FILE *fp = fopen(path, "rb");
    assert(fp);
    char magic[32], actual[2046];
    assert(fgets(magic, sizeof(magic), fp));
    assert(strcmp(magic, "CTOOLS_LABELS_1\n") == 0);
    for (size_t i = 0; i < 5; i++)
    {
        int code;
        size_t length;
        assert(fscanf(fp, "%d %zu", &code, &length) == 2);
        assert(code == codes[i] && length == strlen(texts[i]));
        assert(fgetc(fp) == '\n');
        assert(fread(actual, 1, length, fp) == length);
        assert(memcmp(actual, texts[i], length) == 0);
    }
    assert(fgetc(fp) == EOF && !ferror(fp));
    assert(fclose(fp) == 0);
    unlink(path);
    assert(ctools_label_write_file(texts, codes, 5, "/nonexistent/ctools-labels") == -1);
#ifdef __linux__
    assert(ctools_label_write_file(texts, codes, 5, "/dev/full") == -1);
#endif
    puts("PASS label file framing and failure handling");
    return 0;
}

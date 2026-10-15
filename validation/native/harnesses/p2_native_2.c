
#include <assert.h>
#include "cimport/cimport_encoding.c"
int main(void) {
    const char le[] = "\xff\xfe\0\0a\0\0\0";
    const char be[] = "\0\0\xfe\xff\0\0\0a";
    CImportEncodingDetection a = cimport_detect_encoding(le, sizeof(le)-1);
    CImportEncodingDetection b = cimport_detect_encoding(be, sizeof(be)-1);
    assert(a.encoding == CIMPORT_ENC_UTF32LE && a.bom_length == 4);
    assert(b.encoding == CIMPORT_ENC_UTF32BE && b.bom_length == 4);
    return 0;
}

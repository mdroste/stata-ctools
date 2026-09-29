/* Stateless ReadStat charset conversion using Windows' system code pages.
   Unix builds use system iconv; no third-party DLL is required on Windows. */
#ifndef CTOOLS_ICONV_WIN_H
#define CTOOLS_ICONV_WIN_H
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
typedef struct {
    UINT source, destination;
} *iconv_t;
static inline UINT cio_codepage(const char *name) {
    char s[80];
    size_t j = 0;
    for (size_t i = 0; name[i] && j < sizeof(s) - 1; i++)
        if (isalnum((unsigned char)name[i]))
            s[j++] = (char)toupper((unsigned char)name[i]);
    s[j] = 0;
    if (!strcmp(s, "UTF8"))
        return CP_UTF8;
    if (!strcmp(s, "UTF16") || !strcmp(s, "UTF16LE"))
        return 1200;
    if (!strcmp(s, "UTF16BE"))
        return 1201;
    if (!strcmp(s, "ASCII") || !strcmp(s, "USASCII"))
        return 20127;
    if (!strcmp(s, "MACROMAN") || !strcmp(s, "MACINTOSH"))
        return 10000;
    if (!strcmp(s, "SHIFTJIS") || !strcmp(s, "SJIS"))
        return 932;
    if (!strcmp(s, "BIG5"))
        return 950;
    if (!strcmp(s, "GBK") || !strcmp(s, "GB2312"))
        return 936;
    if (!strcmp(s, "GB18030"))
        return 54936;
    if (!strcmp(s, "EUCJP"))
        return 20932;
    if (!strcmp(s, "EUCKR"))
        return 51949;
    if (!strcmp(s, "ISO2022JP"))
        return 50220;
    if (!strcmp(s, "ISO2022KR"))
        return 50225;
    if (!strcmp(s, "ISO2022CN"))
        return 50227;
    if (!strcmp(s, "KOI8R"))
        return 20866;
    if (!strncmp(s, "ISO8859", 7))
        return 28590 + (UINT)atoi(s + 7);
    if (!strncmp(s, "WINDOWS", 7))
        return (UINT)atoi(s + 7);
    if (!strncmp(s, "CP", 2))
        return (UINT)atoi(s + 2);
    if (!strncmp(s, "IBM", 3))
        return (UINT)atoi(s + 3);
    return 0;
}
static inline iconv_t cio_iconv_open(const char *to, const char *from) {
    UINT source = cio_codepage(from), dest = cio_codepage(to);
    if (!source || !dest || (source != 1200 && source != 1201 && !IsValidCodePage(source)) ||
        (dest != 1200 && dest != 1201 && !IsValidCodePage(dest))) {
        errno = EINVAL;
        return (iconv_t)-1;
    }
    iconv_t c = malloc(sizeof(*c));
    if (!c) {
        errno = ENOMEM;
        return (iconv_t)-1;
    }
    c->source = source;
    c->destination = dest;
    return c;
}
static inline int cio_iconv_close(iconv_t c) {
    free(c);
    return 0;
}
static inline size_t cio_iconv(iconv_t c, char **in, size_t *inleft, char **out, size_t *outleft) {
    if (!in || !*in)
        return 0;
    if (*inleft > INT_MAX) {
        errno = E2BIG;
        return (size_t)-1;
    }
    int n;
    if (c->source == 1200 || c->source == 1201) {
        if (*inleft % 2) {
            errno = EINVAL;
            return (size_t)-1;
        }
        n = (int)(*inleft / 2);
    } else
        n = MultiByteToWideChar(c->source, 0, *in, (int)*inleft, NULL, 0);
    if (!n && *inleft) {
        errno = EILSEQ;
        return (size_t)-1;
    }
    WCHAR *wide = malloc(((size_t)n + 1) * sizeof(WCHAR));
    if (!wide) {
        errno = ENOMEM;
        return (size_t)-1;
    }
    if (c->source == 1200 || c->source == 1201) {
        for (int i = 0; i < n; i++) {
            unsigned char *b = (unsigned char *)*in + 2 * i;
            wide[i] = c->source == 1200 ? (WCHAR)(b[0] | b[1] << 8) : (WCHAR)(b[1] | b[0] << 8);
        }
    } else
        MultiByteToWideChar(c->source, 0, *in, (int)*inleft, wide, n);
    int bytes = c->destination == 1200 || c->destination == 1201
                    ? n * 2
                    : WideCharToMultiByte(c->destination, 0, wide, n, NULL, 0, NULL, NULL);
    if (bytes < 0 || (size_t)bytes > *outleft) {
        free(wide);
        errno = E2BIG;
        return (size_t)-1;
    }
    if (c->destination == 1200 || c->destination == 1201) {
        for (int i = 0; i < n; i++) {
            unsigned char *b = (unsigned char *)*out + 2 * i;
            b[c->destination == 1200 ? 0 : 1] = (unsigned char)wide[i];
            b[c->destination == 1200 ? 1 : 0] = (unsigned char)(wide[i] >> 8);
        }
    } else if (!WideCharToMultiByte(c->destination, 0, wide, n, *out, bytes, NULL, NULL) && n) {
        free(wide);
        errno = EILSEQ;
        return (size_t)-1;
    }
    free(wide);
    *in += *inleft;
    *inleft = 0;
    *out += bytes;
    *outleft -= (size_t)bytes;
    return 0;
}
#define iconv_open cio_iconv_open
#define iconv_close cio_iconv_close
#define iconv cio_iconv
#endif

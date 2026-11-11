/* Stateless ReadStat charset conversion using Windows' system code pages.
   Unix builds use system iconv; no third-party DLL is required on Windows. */
#ifndef CTOOLS_ICONV_WIN_H
#define CTOOLS_ICONV_WIN_H
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
typedef struct
{
    UINT source, destination;
} *iconv_t;
/* An IBM or CP number as Windows numbers the code page (the EBCDIC pages
   Windows has in the 20000s). */
static inline UINT cio_ibm_codepage(UINT n)
{
    static const UINT ebcdic[][2] = {{273, 20273}, {277, 20277}, {278, 20278}, {280, 20280},
                                     {284, 20284}, {285, 20285}, {290, 20290}, {297, 20297},
                                     {420, 20420}, {423, 20423}, {424, 20424}, {838, 20838},
                                     {871, 20871}, {880, 20880}, {905, 20905}, {924, 20924},
                                     {1025, 21025}};
    for (size_t i = 0; i < sizeof(ebcdic) / sizeof(*ebcdic); i++)
        if (ebcdic[i][0] == n)
            return ebcdic[i][1];
    return n;
}
static inline UINT cio_codepage(const char *name)
{
    /* Names compare ignoring case and every character but letters and
       digits (ISO-8859-1, iso8859_1 and ISO88591 are one name). */
    static const struct
    {
        const char *name;
        UINT page;
    } names[] = {{"UTF8", CP_UTF8}, {"UTF16", 1200}, {"UTF16LE", 1200}, {"UTF16BE", 1201},
                 {"ASCII", 20127}, {"USASCII", 20127},
                 {"LATIN1", 28591}, {"L1", 28591}, {"LATIN2", 28592}, {"L2", 28592},
                 {"LATIN3", 28593}, {"L3", 28593}, {"LATIN4", 28594}, {"L4", 28594},
                 {"LATIN5", 28599}, {"L5", 28599}, {"LATIN9", 28605}, {"L9", 28605},
                 {"LATIN0", 28605}, {"CYRILLIC", 28595}, {"ARABIC", 28596}, {"GREEK", 28597},
                 {"HEBREW", 28598}, {"MACROMAN", 10000}, {"MACINTOSH", 10000},
                 {"MACARABIC", 10004}, {"MACHEBREW", 10005}, {"MACGREEK", 10006},
                 {"MACCYRILLIC", 10007}, {"MACROMANIA", 10010}, {"MACUKRAINE", 10017},
                 {"MACTHAI", 10021}, {"MACCENTRALEUROPE", 10029}, {"MACICELAND", 10079},
                 {"MACTURKISH", 10081}, {"MACCROATIAN", 10082}, {"SHIFTJIS", 932}, {"SJIS", 932},
                 {"BIG5", 950}, {"GBK", 936}, {"GB2312", 936}, {"GB18030", 54936},
                 {"EUCJP", 20932}, {"EUCKR", 51949}, {"ISO2022JP", 50220}, {"ISO2022KR", 50225},
                 {"ISO2022CN", 50227}, {"KOI8R", 20866}, {"KOI8U", 21866}, {"TIS620", 874},
                 {"JOHAB", 1361}, {"ISO885911", 874}};
    char s[80];
    size_t j = 0;
    for (size_t i = 0; name[i] && j < sizeof(s) - 1; i++)
        if (isalnum((unsigned char)name[i]))
            s[j++] = (char)toupper((unsigned char)name[i]);
    s[j] = 0;
    for (size_t i = 0; i < sizeof(names) / sizeof(*names); i++)
        if (!strcmp(s, names[i].name))
            return names[i].page;
    if (!strncmp(s, "ISO8859", 7))
        return 28590 + (UINT)atoi(s + 7);
    if (!strncmp(s, "WINDOWS", 7))
        return (UINT)atoi(s + 7);
    if (!strncmp(s, "CP", 2))
        return cio_ibm_codepage((UINT)atoi(s + 2));
    if (!strncmp(s, "IBM", 3))
        return cio_ibm_codepage((UINT)atoi(s + 3));
    return 0;
}
static inline iconv_t cio_iconv_open(const char *to, const char *from)
{
    UINT source = cio_codepage(from), dest = cio_codepage(to);
    if (!source || !dest || (source != 1200 && source != 1201 && !IsValidCodePage(source)) ||
        (dest != 1200 && dest != 1201 && !IsValidCodePage(dest)))
    {
        errno = EINVAL;
        return (iconv_t)-1;
    }
    iconv_t c = malloc(sizeof(*c));
    if (!c)
    {
        errno = ENOMEM;
        return (iconv_t)-1;
    }
    c->source = source;
    c->destination = dest;
    return c;
}
static inline int cio_iconv_close(iconv_t c)
{
    free(c);
    return 0;
}
static inline size_t cio_iconv(iconv_t c, char **in, size_t *inleft, char **out, size_t *outleft)
{
    if (!in || !*in)
        return 0;
    if (*inleft > INT_MAX)
    {
        errno = E2BIG;
        return (size_t)-1;
    }
    /* Like iconv, convert every complete unit and leave an incomplete
       trailing UTF-16 byte in *inleft, reported as EINVAL. */
    size_t usable = *inleft;
    int partial = 0, n;
    if (c->source == 1200 || c->source == 1201)
    {
        partial = (int)(usable % 2);
        usable -= (size_t)partial;
        n = (int)(usable / 2);
    }
    else
        n = MultiByteToWideChar(c->source, 0, *in, (int)usable, NULL, 0);
    if (!n && usable)
    {
        errno = EILSEQ;
        return (size_t)-1;
    }
    WCHAR *wide = malloc(((size_t)n + 1) * sizeof(WCHAR));
    if (!wide)
    {
        errno = ENOMEM;
        return (size_t)-1;
    }
    if (c->source == 1200 || c->source == 1201)
    {
        for (int i = 0; i < n; i++)
        {
            unsigned char *b = (unsigned char *)*in + 2 * i;
            wide[i] = c->source == 1200 ? (WCHAR)(b[0] | b[1] << 8) : (WCHAR)(b[1] | b[0] << 8);
        }
    }
    else
        MultiByteToWideChar(c->source, 0, *in, (int)usable, wide, n);
    int bytes = c->destination == 1200 || c->destination == 1201
                    ? n * 2
                    : WideCharToMultiByte(c->destination, 0, wide, n, NULL, 0, NULL, NULL);
    /* A failed size query (e.g. output over INT_MAX) is E2BIG, not EILSEQ,
       so a caller that skips bad bytes does not retry the whole buffer. */
    if (bytes < 0 || (!bytes && n) || (size_t)bytes > *outleft)
    {
        free(wide);
        errno = E2BIG;
        return (size_t)-1;
    }
    if (c->destination == 1200 || c->destination == 1201)
    {
        for (int i = 0; i < n; i++)
        {
            unsigned char *b = (unsigned char *)*out + 2 * i;
            b[c->destination == 1200 ? 0 : 1] = (unsigned char)wide[i];
            b[c->destination == 1200 ? 1 : 0] = (unsigned char)(wide[i] >> 8);
        }
    }
    else if (!WideCharToMultiByte(c->destination, 0, wide, n, *out, bytes, NULL, NULL) && n)
    {
        free(wide);
        errno = EILSEQ;
        return (size_t)-1;
    }
    free(wide);
    *in += usable;
    *inleft -= usable;
    *out += bytes;
    *outleft -= (size_t)bytes;
    if (partial)
    {
        errno = EINVAL;
        return (size_t)-1;
    }
    return 0;
}
#define iconv_open cio_iconv_open
#define iconv_close cio_iconv_close
#define iconv cio_iconv
#endif

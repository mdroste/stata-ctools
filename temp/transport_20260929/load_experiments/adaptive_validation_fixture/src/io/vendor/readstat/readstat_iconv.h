#ifndef READSTAT_ICONV_H
#define READSTAT_ICONV_H
#ifdef _WIN32
#include "../../cio_iconv_win.h"
#else
#include <iconv.h>
#endif

/* ICONV_CONST defined by autotools during configure according
 * to the current platform. Some people copy-paste the source code, so
 * provide some fallback logic */
#ifndef ICONV_CONST
#define ICONV_CONST
#endif

typedef ICONV_CONST char ** readstat_iconv_inbuf_t;

typedef struct readstat_charset_entry_s {
    int     code;
    char    name[32];
} readstat_charset_entry_t;

#endif

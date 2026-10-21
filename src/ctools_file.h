#ifndef CTOOLS_FILE_H
#define CTOOLS_FILE_H
/* All plugin-owned paths are UTF-8. Never change the host's locale/code page. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <wchar.h>

static inline wchar_t *ctools_wpath(const char *path)
{
    int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
    if (!n) { errno = EINVAL; return NULL; }
    wchar_t *input = malloc((size_t)n * sizeof(wchar_t));
    if (!input) { SetLastError(ERROR_NOT_ENOUGH_MEMORY); return NULL; }
    if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, input, n)) {
        free(input); errno = EINVAL; return NULL;
    }
    for (wchar_t *p = input; *p; p++) if (*p == L'/') *p = L'\\';
    if (!wcsncmp(input, L"\\\\?\\", 4)) return input;
    DWORD need = GetFullPathNameW(input, 0, NULL, NULL);
    wchar_t *absolute = need ? malloc(((size_t)need + 8) * sizeof(wchar_t)) : NULL;
    if (!absolute) { free(input); return NULL; }
    DWORD got = GetFullPathNameW(input, need, absolute + 8, NULL);
    free(input);
    if (!got || got >= need) { free(absolute); errno = EINVAL; return NULL; }
    if (absolute[8] == L'\\' && absolute[9] == L'\\') {
        memmove(absolute + 8, absolute + 10, (wcslen(absolute + 10) + 1) * sizeof(wchar_t));
        memcpy(absolute, L"\\\\?\\UNC\\", 8 * sizeof(wchar_t));
    } else {
        memmove(absolute + 4, absolute + 8, (wcslen(absolute + 8) + 1) * sizeof(wchar_t));
        memcpy(absolute, L"\\\\?\\", 4 * sizeof(wchar_t));
    }
    return absolute;
}
static inline FILE *ctools_fopen(const char *path, const char *mode)
{
    wchar_t *wide = ctools_wpath(path), wm[16];
    if (!wide) return NULL;
    size_t n = strlen(mode);
    if (n >= 16) { free(wide); errno = EINVAL; return NULL; }
    for (size_t i = 0; i <= n; i++) wm[i] = (wchar_t)(unsigned char)mode[i];
    FILE *fp = _wfopen(wide, wm); free(wide); return fp;
}
static inline int ctools_remove(const char *path)
{
    wchar_t *wide = ctools_wpath(path); if (!wide) return -1;
    int rc = _wremove(wide); free(wide); return rc;
}
static inline HANDLE ctools_create_file(const char *path, DWORD access, DWORD share,
                                       DWORD creation, DWORD flags)
{
    wchar_t *wide = ctools_wpath(path); if (!wide) return INVALID_HANDLE_VALUE;
    HANDLE h = CreateFileW(wide, access, share, NULL, creation, flags, NULL);
    DWORD error = GetLastError(); free(wide); SetLastError(error); return h;
}
static inline int ctools_move_file(const char *from, const char *to, int replace)
{
    wchar_t *a = ctools_wpath(from), *b = ctools_wpath(to);
    if (!a || !b) { free(a); free(b); return -1; }
    BOOL ok = MoveFileExW(a, b, MOVEFILE_WRITE_THROUGH | (replace ? MOVEFILE_REPLACE_EXISTING : 0));
    DWORD error = GetLastError(); free(a); free(b); SetLastError(error);
    if (!ok) errno = error == ERROR_ALREADY_EXISTS || error == ERROR_FILE_EXISTS ? EEXIST : EIO;
    return ok ? 0 : -1;
}
#define ctools_ftell _ftelli64
#define ctools_fseek _fseeki64
#else
#include <unistd.h>
#define ctools_fopen fopen
#define ctools_remove remove
#define ctools_ftell ftello
#define ctools_fseek fseeko
#endif
#endif

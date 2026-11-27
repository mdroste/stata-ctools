/*
 * csql_api.c - load the DuckDB shared library and resolve its C API.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include "csql/csql_api.h"

csql_api csql_db;

#ifdef _WIN32
static void *csql_lib_open(const char *path, char *msg, size_t cap)
{
    int n = MultiByteToWideChar(CP_UTF8, 0, path, -1, NULL, 0);
    wchar_t *wide = n > 0 ? (wchar_t *)malloc((size_t)n * sizeof(wchar_t)) : NULL;
    HMODULE h = NULL;
    if (wide && MultiByteToWideChar(CP_UTF8, 0, path, -1, wide, n) > 0) {
        /* The full path, and the library's folder for any dependency. */
        h = LoadLibraryExW(wide, NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
    }
    if (!h) {
        DWORD err = GetLastError();
        char text[256] = "";
        FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, err, 0,
                       text, (DWORD)sizeof(text), NULL);
        size_t len = strlen(text);
        while (len && (text[len - 1] == '\n' || text[len - 1] == '\r' || text[len - 1] == ' '))
            text[--len] = '\0';
        snprintf(msg, cap, "could not load %s: %s (Windows error %lu)", path, text,
                 (unsigned long)err);
    }
    free(wide);
    return (void *)h;
}

static void *csql_lib_symbol(void *lib, const char *name)
{
    return (void *)GetProcAddress((HMODULE)lib, name);
}
#else
static void *csql_lib_open(const char *path, char *msg, size_t cap)
{
    void *h = dlopen(path, RTLD_NOW | RTLD_LOCAL);
    if (!h) {
        const char *err = dlerror();
        snprintf(msg, cap, "could not load %s: %s", path, err ? err : "unknown error");
    }
    return h;
}

static void *csql_lib_symbol(void *lib, const char *name)
{
    return dlsym(lib, name);
}
#endif

int csql_api_load(const char *path, char *msg, size_t cap)
{
    csql_api api;
    void *lib;
    if (msg && cap) msg[0] = '\0';
    if (!path || !*path) {
        snprintf(msg, cap, "no DuckDB library given");
        return 601;
    }
    if (csql_db.loaded && strcmp(csql_db.path, path) == 0) return 0;
    lib = csql_lib_open(path, msg, cap);
    if (!lib) return 601;

    memset(&api, 0, sizeof(api));
    /* Function pointers and void * share a representation on every platform
     * ctools supports (dlsym and GetProcAddress rely on it). */
#define CSQL_API_RESOLVE(name, required)                                                           \
    {                                                                                              \
        void *sym = csql_lib_symbol(lib, #name);                                                   \
        memcpy(&api.name, &sym, sizeof(sym));                                                      \
        if (!sym && (required)) {                                                                  \
            snprintf(msg, cap,                                                                     \
                     "%s is not a usable DuckDB library: it has no function " #name                \
                     " (csql needs DuckDB 1.1 or later)",                                          \
                     path);                                                                        \
            return 610;                                                                            \
        }                                                                                          \
    }
    CSQL_API_FUNCTIONS(CSQL_API_RESOLVE)
#undef CSQL_API_RESOLVE

    /* The library stays loaded: DuckDB registers thread-local and static
     * state that must outlive any database closed in this session. */
    snprintf(api.path, sizeof(api.path), "%s", path);
    {
        const char *v = api.duckdb_library_version();
        snprintf(api.version, sizeof(api.version), "%s", v ? v : "");
    }
    api.loaded = 1;
    csql_db = api;
    return 0;
}

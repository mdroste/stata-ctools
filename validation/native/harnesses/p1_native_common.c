
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
static int calls, fail_at, live;
static void *fi_malloc(size_t n) {
    if (++calls == fail_at) return NULL;
    void *p = malloc(n);
    if (p) { memset(p, 0xa5, n); live++; }
    return p;
}
static void *fi_calloc(size_t n, size_t size) {
    if (++calls == fail_at) return NULL;
    void *p = calloc(n, size);
    if (p) live++;
    return p;
}
static int fi_aligned(void **p, size_t align, size_t n) {
    if (++calls == fail_at) return ENOMEM;
    int rc = posix_memalign(p, align, n);
    if (!rc) live++;
    return rc;
}
static void fi_free(void *p) { if (p) live--; free(p); }
#define malloc fi_malloc
#define calloc fi_calloc
#define posix_memalign fi_aligned
#define free fi_free

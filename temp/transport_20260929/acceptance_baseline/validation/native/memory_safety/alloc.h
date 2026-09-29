/* Failure injection and ownership tracking for native command regressions. */
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static size_t live_bytes, peak_bytes;
typedef struct {
    void *p;
    size_t n;
    const char *file;
    int line;
} Entry;
static Entry blocks[100000];
static int calls, fail_at, alive;
static void record(void *p, size_t n, const char *f, int l) {
    if (!p)
        return;
    for (int i = 0; i < 100000; i++)
        if (!blocks[i].p) {
            blocks[i] = (Entry){p, n, f, l};
            alive++;
            live_bytes += n;
            if (live_bytes > peak_bytes)
                peak_bytes = live_bytes;
            return;
        }
    abort();
}
static void *tm_impl(size_t n, const char *f, int l) {
    if (++calls == fail_at)
        return NULL;
    void *p = malloc(n);
    if (p)
        memset(p, 0xa5, n);
    record(p, n, f, l);
    return p;
}
static void *tc_impl(size_t n, size_t s, const char *f, int l) {
    if (++calls == fail_at)
        return NULL;
    void *p = calloc(n, s);
    record(p, n * s, f, l);
    return p;
}
static void tf_impl(void *p) {
    if (!p)
        return;
    for (int i = 0; i < 100000; i++)
        if (blocks[i].p == p) {
            live_bytes -= blocks[i].n;
            blocks[i].p = NULL;
            alive--;
            free(p);
            return;
        }
    fprintf(stderr, "UNTRACKED FREE %p\n", p);
    abort();
}
static int ta_impl(void **p, size_t a, size_t n, const char *f, int l) {
    if (++calls == fail_at)
        return ENOMEM;
    int rc = posix_memalign(p, a, n);
    if (!rc)
        record(*p, n, f, l);
    return rc;
}
static void leftovers(void) {
    for (int i = 0; i < 100000; i++)
        if (blocks[i].p)
            printf("LEAK %zu %s:%d\n", blocks[i].n, blocks[i].file, blocks[i].line);
}
static void *tr_impl(void *p, size_t n, const char *f, int l) {
    if (++calls == fail_at)
        return NULL;
    int idx = -1;
    for (int i = 0; i < 100000; i++)
        if (p && blocks[i].p == p) {
            idx = i;
            break;
        }
    void *q = realloc(p, n);
    if (q) {
        if (idx >= 0) {
            live_bytes -= blocks[idx].n;
            blocks[idx].p = NULL;
            alive--;
        }
        record(q, n, f, l);
    }
    return q;
}

#include <pthread.h>
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
static void *tm(size_t n, const char *f, int l) {
    pthread_mutex_lock(&alloc_mutex);
    void *result = tm_impl(n, f, l);
    pthread_mutex_unlock(&alloc_mutex);
    return result;
}
static void *tc(size_t n, size_t z, const char *f, int l) {
    pthread_mutex_lock(&alloc_mutex);
    void *result = tc_impl(n, z, f, l);
    pthread_mutex_unlock(&alloc_mutex);
    return result;
}
static void *tr(void *p, size_t n, const char *f, int l) {
    pthread_mutex_lock(&alloc_mutex);
    void *result = tr_impl(p, n, f, l);
    pthread_mutex_unlock(&alloc_mutex);
    return result;
}
static int ta(void **p, size_t a, size_t n, const char *f, int l) {
    pthread_mutex_lock(&alloc_mutex);
    int result = ta_impl(p, a, n, f, l);
    pthread_mutex_unlock(&alloc_mutex);
    return result;
}
static void tf(void *p) {
    pthread_mutex_lock(&alloc_mutex);
    tf_impl(p);
    pthread_mutex_unlock(&alloc_mutex);
}
#define realloc(p, n) tr(p, n, __FILE__, __LINE__)
#define malloc(n) tm(n, __FILE__, __LINE__)
#define calloc(n, s) tc(n, s, __FILE__, __LINE__)
#define free(p) tf(p)
#define posix_memalign(p, a, n) ta(p, a, n, __FILE__, __LINE__)

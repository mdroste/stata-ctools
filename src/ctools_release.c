/*
 * ctools_release.c - return a finished command's memory on a helper thread
 * (see ctools_release.h)
 *
 * Only the thread that runs stata_call() opens a session, queues, starts and
 * joins the helper; a worker thread that reaches a release call frees at
 * once. The image is pinned (ctools_plugin.c), so the helper's code stays
 * mapped even if Stata unloads the plugin while it runs.
 */

#include "ctools_release.h"

#include <pthread.h>
#include <string.h>

typedef struct {
    void *p;
    ctools_release_fn fn;
    int kind;
} release_entry;

typedef struct {
    release_entry *items;
    size_t n, cap;
} release_queue;

static _Thread_local int session_open;
static release_queue pending;   /* filled during a session */
static release_queue handed;    /* owned by the helper while it runs */
static pthread_t helper;
static int helper_running;

int ctools_release_active(void)
{
    return session_open;
}

int ctools_release_queue(void *p, int kind, ctools_release_fn fn)
{
    if (!session_open) return 0;
    if (pending.n == pending.cap) {
        size_t cap = pending.cap ? 2 * pending.cap : 64;
        release_entry *grown = (release_entry *)realloc(pending.items, cap * sizeof(*grown));
        if (!grown) return 0;
        pending.items = grown;
        pending.cap = cap;
    }
    pending.items[pending.n].p = p;
    pending.items[pending.n].fn = fn;
    pending.items[pending.n].kind = kind;
    pending.n++;
    return 1;
}

static void release_run(release_queue *q)
{
    for (size_t i = 0; i < q->n; i++) {
        release_entry *e = &q->items[i];
        if (e->kind == CTOOLS_RELEASE_CALL) e->fn(e->p);
        else if (e->kind == CTOOLS_RELEASE_ALIGNED) ctools_aligned_free(e->p);
        else free(e->p);
    }
    free(q->items);
    memset(q, 0, sizeof(*q));
}

static void (*helper_then)(void);  /* set before the helper starts */

static void *release_main(void *arg)
{
    release_run((release_queue *)arg);
    if (helper_then) helper_then();
    return NULL;
}

void ctools_release_wait(void)
{
    if (helper_running) {
        pthread_join(helper, NULL);
        helper_running = 0;
    }
}

void ctools_release_open(void)
{
    session_open = 1;
}

void ctools_release_end(void)
{
    session_open = 0;
}

void ctools_release_close(void (*then)(void))
{
    session_open = 0;
    if (pending.n == 0) {
        if (then) then();
        return;
    }
    ctools_release_wait();
    handed = pending;
    memset(&pending, 0, sizeof(pending));
    helper_then = then;
    if (pthread_create(&helper, NULL, release_main, &handed) == 0) {
        helper_running = 1;
        return;
    }
    release_run(&handed);
    if (then) then();
}

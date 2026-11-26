/*
 * ctools_release.h - return a finished command's memory on a helper thread
 *
 * Releasing the large blocks of a finished estimation takes tens of
 * milliseconds on Windows (40-90 ms for an HDFE fit on 10 million rows),
 * time that Stata used to wait for before it could post the results. During
 * a plugin call, stata_call() opens a release session: a command's final
 * cleanup hands its blocks to ctools_release_later*(), which queue them, and
 * when the command returns a helper thread releases the queue while Stata
 * posts and displays the results. Every later plugin call except the
 * version and clock queries waits for the helper before it starts, so the
 * memory of one command is returned before the next one allocates.
 *
 * Outside a session (native harnesses call the modules directly), and when
 * the queue cannot grow, each call releases the block at once, where it was
 * released before. The calls are macros so that this immediate release uses
 * the deallocator in effect where they are written (harnesses redefine it).
 *
 * Use them only in a command's final cleanup, which runs once per plugin
 * call: a queued block stays allocated until the call returns.
 */

#ifndef CTOOLS_RELEASE_H
#define CTOOLS_RELEASE_H

#include <stdlib.h>
#include "ctools_alloc.h"

typedef void (*ctools_release_fn)(void *arg);

enum { CTOOLS_RELEASE_FREE, CTOOLS_RELEASE_ALIGNED, CTOOLS_RELEASE_CALL };

/* 1 when the entry was queued, 0 when the caller must release it now */
int ctools_release_queue(void *p, int kind, ctools_release_fn fn);
/* 1 inside a release session, on the thread that opened it */
int ctools_release_active(void);

/* free(p), on the helper or now */
#define ctools_release_later(p) do { \
        void *ctools_release_p_ = (void *)(p); \
        if (ctools_release_p_ && \
            !ctools_release_queue(ctools_release_p_, CTOOLS_RELEASE_FREE, NULL)) \
            free(ctools_release_p_); \
    } while (0)

/* ctools_aligned_free(p), on the helper or now */
#define ctools_release_later_aligned(p) do { \
        void *ctools_release_p_ = (void *)(p); \
        if (ctools_release_p_ && \
            !ctools_release_queue(ctools_release_p_, CTOOLS_RELEASE_ALIGNED, NULL)) \
            ctools_aligned_free(ctools_release_p_); \
    } while (0)

/* fn(arg), on the helper or now: a cleanup that only frees memory that arg
 * owns, and that nothing else refers to any more */
#define ctools_release_later_call(fn, arg) do { \
        void *ctools_release_a_ = (void *)(arg); \
        if (!ctools_release_queue(ctools_release_a_, CTOOLS_RELEASE_CALL, (fn))) \
            (fn)(ctools_release_a_); \
    } while (0)

/* For stata_call() only: wait for the helper before a command; open the
 * session for the command's handler and end it when the handler returns;
 * close when the call is done: a helper releases the queue and then runs
 * then(), or then() runs at once when nothing was queued */
void ctools_release_wait(void);
void ctools_release_open(void);
void ctools_release_end(void);
void ctools_release_close(void (*then)(void));

#endif /* CTOOLS_RELEASE_H */

#include "ctools_cancel.h"
#ifdef _OPENMP
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif
#endif

static atomic_int cancelled;
static _Thread_local int calling_thread;
static ctools_cancel_poll host_poll;

void ctools_cancel_begin(ctools_cancel_poll poll)
{
    host_poll = poll;
    calling_thread = 1;
    atomic_store_explicit(&cancelled, 0, memory_order_relaxed);
}

void ctools_cancel_end(void)
{
    calling_thread = 0;
    host_poll = 0;
    atomic_store_explicit(&cancelled, 0, memory_order_relaxed);
}

int ctools_cancel_requested(void)
{
    if (atomic_load_explicit(&cancelled, memory_order_relaxed))
        return 1;
    if (calling_thread && host_poll && host_poll())
    {
        atomic_store_explicit(&cancelled, 1, memory_order_relaxed);
        return 1;
    }
    return 0;
}

void ctools_cancel_wait(ctools_cancel_team *team)
{
#ifdef _OPENMP
    int size = omp_get_num_threads();
    atomic_fetch_add_explicit(&team->arrived, 1, memory_order_release);
    while (atomic_load_explicit(&team->arrived, memory_order_acquire) < size)
    {
        ctools_cancel_requested();
#ifdef _WIN32
        Sleep(1);
#else
        const struct timespec delay = {0, 1000000};
        nanosleep(&delay, 0);
#endif
    }
#else
    (void)team;
#endif
    ctools_cancel_requested();
}

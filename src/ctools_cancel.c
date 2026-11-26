#include "ctools_cancel.h"
#ifdef _OPENMP
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#include <time.h>
#endif
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <emmintrin.h>
#define cancel_pause() _mm_pause()
#elif defined(__aarch64__) || defined(_M_ARM64)
#define cancel_pause() __asm__ __volatile__("yield")
#else
#define cancel_pause() ((void)0)
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

#ifdef _OPENMP
/*
 * Where the calling thread sleeps while its team works. One lock and one
 * condition serve every team: only the command's calling thread ever sleeps
 * here, and a wake-up meant for another team's counter only makes it look at
 * its own again.
 *
 * It used to sleep a millisecond at a time, and so did every other member
 * that finished early, although only the calling thread can poll the host.
 * A parallel region then ended on the slowest sleeper: about a millisecond
 * late where sleeps are that fine, and up to a timer tick late on Windows
 * (15.6 ms unless the process asked for a finer one). A projection on one
 * fixed effect of a thousand rows took 15 ms with sixteen threads and
 * microseconds with one.
 */
#ifdef _WIN32
static SRWLOCK team_lock_ = SRWLOCK_INIT;
static CONDITION_VARIABLE team_done_ = CONDITION_VARIABLE_INIT;
static void team_lock(void) { AcquireSRWLockExclusive(&team_lock_); }
static void team_unlock(void) { ReleaseSRWLockExclusive(&team_lock_); }
static void team_wake(void) { WakeAllConditionVariable(&team_done_); }
/* With the lock held: sleep until woken, or for about a millisecond */
static void team_sleep(void) { SleepConditionVariableSRW(&team_done_, &team_lock_, 1, 0); }
#else
static pthread_mutex_t team_lock_ = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t team_done_ = PTHREAD_COND_INITIALIZER;
static void team_lock(void) { pthread_mutex_lock(&team_lock_); }
static void team_unlock(void) { pthread_mutex_unlock(&team_lock_); }
static void team_wake(void) { pthread_cond_broadcast(&team_done_); }
static void team_sleep(void)
{
    struct timespec until;
    clock_gettime(CLOCK_REALTIME, &until);
    until.tv_nsec += 1000000L;
    if (until.tv_nsec >= 1000000000L)
    {
        until.tv_sec += 1;
        until.tv_nsec -= 1000000000L;
    }
    pthread_cond_timedwait(&team_done_, &team_lock_, &until);
}
#endif

/* Looks at the counter the calling thread takes before it sleeps: a team
 * that finishes together costs it neither a sleep nor a wake-up */
#define CANCEL_SPINS 256
#endif

void ctools_cancel_wait(ctools_cancel_team *team)
{
#ifdef _OPENMP
    const int size = omp_get_num_threads();
    const int arrived = atomic_fetch_add(&team->arrived, 1) + 1;
    if (!calling_thread)
    {
        /* Only the calling thread can poll the host: any other member goes
         * on to the region's barrier. The last one in wakes the calling
         * thread, which may be asleep waiting for this team. Taking the lock
         * puts the wake-up either before the sleeper's look at the counter
         * or after it has gone to sleep, never in between. */
        if (arrived >= size)
        {
            team_lock();
            team_wake();
            team_unlock();
        }
        return;
    }
    for (int spin = 0; spin < CANCEL_SPINS && atomic_load(&team->arrived) < size; spin++)
        cancel_pause();
    for (;;)
    {
        team_lock();
        const int complete = atomic_load(&team->arrived) >= size;
        if (!complete)
            team_sleep();
        team_unlock();
        if (complete || atomic_load(&team->arrived) >= size)
            break;
        /* still working: poll the host, outside the lock */
        ctools_cancel_requested();
    }
#else
    (void)team;
#endif
    ctools_cancel_requested();
}

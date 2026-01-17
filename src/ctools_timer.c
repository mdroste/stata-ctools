/*
 * ctools_timer.c - Cross-platform high-resolution timing utilities
 *
 * Platform detection and implementation:
 * - macOS: mach_absolute_time (nanosecond resolution)
 * - Linux: clock_gettime with CLOCK_MONOTONIC (nanosecond resolution)
 * - Windows: QueryPerformanceCounter (sub-microsecond resolution)
 * - Fallback: gettimeofday (microsecond resolution)
 *
 * All implementations provide at least millisecond accuracy.
 */

#include "ctools_timer.h"

/* Platform detection */
#if defined(__APPLE__) && defined(__MACH__)
    #define CTOOLS_TIMER_MACH 1
#elif defined(_WIN32) || defined(_WIN64)
    #define CTOOLS_TIMER_WINDOWS 1
#elif defined(__linux__) || defined(__unix__) || defined(_POSIX_VERSION)
    #define CTOOLS_TIMER_POSIX 1
#else
    #define CTOOLS_TIMER_FALLBACK 1
#endif

/* Platform-specific includes and state */
#if defined(CTOOLS_TIMER_MACH)
    #include <mach/mach_time.h>
    #include <pthread.h>
    static mach_timebase_info_data_t g_timebase_info;
    static pthread_once_t g_timer_once = PTHREAD_ONCE_INIT;

#elif defined(CTOOLS_TIMER_WINDOWS)
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    static LARGE_INTEGER g_frequency;
    static INIT_ONCE g_timer_once = INIT_ONCE_STATIC_INIT;

#elif defined(CTOOLS_TIMER_POSIX)
    #include <time.h>
    #include <pthread.h>
    static pthread_once_t g_timer_once = PTHREAD_ONCE_INIT;

#else /* Fallback */
    #include <sys/time.h>
    #include <pthread.h>
    static pthread_once_t g_timer_once = PTHREAD_ONCE_INIT;
#endif

#if !defined(CTOOLS_TIMER_WINDOWS)
/* Internal initialization function (called once via pthread_once) */
static void ctools_timer_init_internal(void) {
#if defined(CTOOLS_TIMER_MACH)
    mach_timebase_info(&g_timebase_info);
#else
    /* No initialization needed for clock_gettime or gettimeofday */
#endif
}
#endif /* !CTOOLS_TIMER_WINDOWS */

#if defined(CTOOLS_TIMER_WINDOWS)
/* Windows callback for InitOnceExecuteOnce */
static BOOL CALLBACK ctools_timer_init_callback(
    PINIT_ONCE InitOnce, PVOID Parameter, PVOID *Context) {
    (void)InitOnce; (void)Parameter; (void)Context;
    QueryPerformanceFrequency(&g_frequency);
    return TRUE;
}
#endif

/*
 * Initialize the timer subsystem.
 * Thread-safe: uses pthread_once (POSIX) or InitOnceExecuteOnce (Windows).
 */
void ctools_timer_init(void) {
#if defined(CTOOLS_TIMER_WINDOWS)
    InitOnceExecuteOnce(&g_timer_once, ctools_timer_init_callback, NULL, NULL);
#else
    pthread_once(&g_timer_once, ctools_timer_init_internal);
#endif
}

/*
 * Get current time in seconds.
 */
double ctools_timer_seconds(void) {
    /* Thread-safe auto-initialize */
    ctools_timer_init();

#if defined(CTOOLS_TIMER_MACH)
    uint64_t t = mach_absolute_time();
    /* Convert to nanoseconds, then to seconds */
    double nanos = (double)t * g_timebase_info.numer / g_timebase_info.denom;
    return nanos / 1e9;

#elif defined(CTOOLS_TIMER_WINDOWS)
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)g_frequency.QuadPart;

#elif defined(CTOOLS_TIMER_POSIX)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;

#else /* Fallback: gettimeofday */
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
#endif
}

/*
 * Get current time in milliseconds.
 */
double ctools_timer_ms(void) {
    return ctools_timer_seconds() * 1000.0;
}

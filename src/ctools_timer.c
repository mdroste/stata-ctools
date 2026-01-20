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
 *
 * Note: Uses simple volatile flags for initialization to avoid
 * static initializer issues with cross-compiled Windows DLLs.
 */

#include "ctools_timer.h"

/* Platform detection */
#if defined(__APPLE__) && defined(__MACH__)
    #define CTOOLS_TIMER_MACH 1
    #include <mach/mach_time.h>
#elif defined(_WIN32) || defined(_WIN64)
    #define CTOOLS_TIMER_WINDOWS 1
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#elif defined(__linux__) || defined(__unix__) || defined(_POSIX_VERSION)
    #define CTOOLS_TIMER_POSIX 1
    #include <time.h>
#else
    #define CTOOLS_TIMER_FALLBACK 1
    #include <sys/time.h>
#endif

/* Platform-specific state - use simple volatile flags, no fancy initializers */
#if defined(CTOOLS_TIMER_MACH)
    static mach_timebase_info_data_t g_timebase_info;
    static volatile int g_timer_initialized = 0;

#elif defined(CTOOLS_TIMER_WINDOWS)
    static LARGE_INTEGER g_frequency;
    static volatile int g_timer_initialized = 0;

#else
    /* POSIX and fallback don't need initialization */
    static volatile int g_timer_initialized = 1;
#endif

/*
 * Initialize the timer subsystem.
 * Simple double-checked pattern - safe enough for our use case since
 * worst case is redundant initialization, not corruption.
 */
void ctools_timer_init(void) {
    if (g_timer_initialized) {
        return;
    }

#if defined(CTOOLS_TIMER_MACH)
    mach_timebase_info(&g_timebase_info);
    g_timer_initialized = 1;

#elif defined(CTOOLS_TIMER_WINDOWS)
    QueryPerformanceFrequency(&g_frequency);
    g_timer_initialized = 1;

#else
    /* Nothing to initialize */
    g_timer_initialized = 1;
#endif
}

/*
 * Get current time in seconds.
 */
double ctools_timer_seconds(void) {
    /* Auto-initialize on first call */
    if (!g_timer_initialized) {
        ctools_timer_init();
    }

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

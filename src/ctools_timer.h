/*
 * ctools_timer.h - Cross-platform high-resolution timing utilities
 *
 * Provides millisecond-accurate timing functions that work on:
 * - macOS (mach_absolute_time)
 * - Linux (clock_gettime with CLOCK_MONOTONIC)
 * - Windows (QueryPerformanceCounter)
 * - Fallback (gettimeofday)
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_TIMER_H
#define CTOOLS_TIMER_H

#include <stdint.h>

/*
 * Initialize the timer subsystem.
 * Call once at program start. Safe to call multiple times.
 */
void ctools_timer_init(void);

/*
 * Get current time in seconds (double precision).
 * Resolution: ~1 microsecond on most platforms.
 * Returns: seconds since an arbitrary epoch (use for elapsed time only)
 */
double ctools_timer_seconds(void);

/*
 * Get current time in milliseconds (double precision).
 * Resolution: ~1 microsecond on most platforms.
 * Returns: milliseconds since an arbitrary epoch (use for elapsed time only)
 */
double ctools_timer_ms(void);

/*
 * Convenience macro for timing a code block.
 * Usage:
 *   double elapsed;
 *   CTOOLS_TIME_BLOCK(elapsed) {
 *       // code to time
 *   }
 *   printf("Elapsed: %.3f ms\n", elapsed);
 */
#define CTOOLS_TIME_BLOCK(elapsed_ms) \
    for (double _t_start = ctools_timer_ms(), _t_done = 0; \
         !_t_done; \
         elapsed_ms = ctools_timer_ms() - _t_start, _t_done = 1)

#endif /* CTOOLS_TIMER_H */

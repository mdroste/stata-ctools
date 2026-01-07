/*
    ctools_config.h
    Configuration parameters for ctools performance tuning

    These values can be adjusted based on profiling results and target hardware.
*/

#ifndef CTOOLS_CONFIG_H
#define CTOOLS_CONFIG_H

/* ============================================================================
   Parallelization thresholds
   ============================================================================ */

/*
    Minimum total bytes to justify parallel I/O (load/store).
    Below this threshold, thread creation overhead exceeds benefit.
    Default: 1MB (assumes ~100us thread creation overhead)
*/
#define MIN_PARALLEL_BYTES (1 * 1024 * 1024)

/*
    Minimum observations per thread for parallel sort.
    Below this, sequential radix sort is faster due to better cache locality.
*/
#define MIN_OBS_PER_THREAD 100000

/*
    Maximum number of threads for parallel operations.
    Usually matches CPU core count; can be tuned for hyperthreading.
*/
#define NUM_THREADS 8

/* ============================================================================
   Memory alignment
   ============================================================================ */

/*
    Cache line size for aligned allocations.
    64 bytes is standard for x86-64; ARM may differ.
*/
#define CACHE_LINE_SIZE 64

/* ============================================================================
   Radix sort parameters
   ============================================================================ */

/*
    Radix size: 8 bits = 256 buckets.
    Larger radix = fewer passes but more memory; 8 is optimal for most cases.
*/
#define RADIX_BITS 8
#define RADIX_SIZE (1 << RADIX_BITS)
#define RADIX_MASK (RADIX_SIZE - 1)

/* ============================================================================
   Prefetch tuning
   ============================================================================ */

/*
    Prefetch distance for indirect memory access patterns.
    Tune based on memory latency and cache hierarchy.
*/
#define PREFETCH_DISTANCE 16

#endif /* CTOOLS_CONFIG_H */

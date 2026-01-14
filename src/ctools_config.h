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

/*
    Unified prefetch macros for all ctools modules.
    CTOOLS_PREFETCH(addr)   - Prefetch for reading
    CTOOLS_PREFETCH_W(addr) - Prefetch for writing
*/
#if defined(__GNUC__) || defined(__clang__)
    #define CTOOLS_PREFETCH(addr)   __builtin_prefetch((addr), 0, 1)
    #define CTOOLS_PREFETCH_W(addr) __builtin_prefetch((addr), 1, 1)
    #define CTOOLS_LIKELY(x)        __builtin_expect(!!(x), 1)
    #define CTOOLS_UNLIKELY(x)      __builtin_expect(!!(x), 0)
#else
    #define CTOOLS_PREFETCH(addr)   ((void)0)
    #define CTOOLS_PREFETCH_W(addr) ((void)0)
    #define CTOOLS_LIKELY(x)        (x)
    #define CTOOLS_UNLIKELY(x)      (x)
#endif

/* Restrict keyword for aliasing hints */
#ifdef __GNUC__
    #define CTOOLS_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define CTOOLS_RESTRICT __restrict
#else
    #define CTOOLS_RESTRICT
#endif

/* Always inline hint */
#if defined(__GNUC__) || defined(__clang__)
    #define CTOOLS_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
    #define CTOOLS_INLINE __forceinline
#else
    #define CTOOLS_INLINE inline
#endif

/* ============================================================================
   Cross-platform aligned memory allocation

   IMPORTANT: On Windows, _aligned_malloc() requires _aligned_free().
   On POSIX systems, posix_memalign() memory can be freed with regular free().
   Always use ctools_aligned_free() for memory allocated with ctools_aligned_alloc().
   ============================================================================ */

#if defined(_WIN32)
    #include <malloc.h>
#else
    #include <stdlib.h>
#endif

/*
    Allocate memory aligned to specified boundary.
    Returns NULL on failure.
    MUST be freed with ctools_aligned_free().
*/
static inline void *ctools_aligned_alloc(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, alignment);
#else
    if (posix_memalign(&ptr, alignment, size) != 0) {
        return NULL;
    }
#endif
    return ptr;
}

/*
    Free memory allocated with ctools_aligned_alloc().
    Safe to call with NULL.
*/
static inline void ctools_aligned_free(void *ptr)
{
    if (ptr == NULL) return;
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ============================================================================
   Cross-platform memory barriers

   IMPORTANT: __sync_synchronize() is a GCC builtin that may not work on
   Windows with MSVC. Use ctools_memory_barrier() for portable code.
   ============================================================================ */

#if defined(_WIN32) && defined(_MSC_VER)
    #include <intrin.h>
    #define ctools_memory_barrier() _ReadWriteBarrier(); MemoryBarrier()
#elif defined(__GNUC__) || defined(__clang__)
    #define ctools_memory_barrier() __sync_synchronize()
#else
    /* Fallback: compiler barrier only */
    #define ctools_memory_barrier() do { __asm__ __volatile__("" ::: "memory"); } while(0)
#endif

#endif /* CTOOLS_CONFIG_H */

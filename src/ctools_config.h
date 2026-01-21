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
#define NUM_THREADS 12

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
   OpenMP compatibility layer
   Provides fallback stubs when OpenMP is not available.
   ============================================================================ */

#ifdef _OPENMP
    #include <omp.h>
    #define CTOOLS_OPENMP_ENABLED 1
#else
    /* Fallback stubs for OpenMP runtime functions */
    static inline int omp_get_thread_num(void)  { return 0; }
    static inline int omp_get_max_threads(void) { return 1; }
    static inline int omp_get_num_threads(void) { return 1; }
    static inline void omp_set_num_threads(int n) { (void)n; }
    #define CTOOLS_OPENMP_ENABLED 0
#endif

/* ============================================================================
   Thread Diagnostics

   Helper macro to save thread information to Stata scalars for verbose output.
   Call at the start of each command to report threading configuration.

   Parameters:
     prefix - scalar name prefix (e.g., "_csort" produces "_csort_threads_max")

   Saves the following scalars:
     <prefix>_threads_max     - Maximum threads available (omp_get_max_threads)
     <prefix>_openmp_enabled  - 1 if OpenMP is enabled, 0 otherwise
   ============================================================================ */

/* Forward declaration - stplugin.h provides SF_scal_save */
#ifndef SF_scal_save
    /* Will be defined when stplugin.h is included */
#endif

#define CTOOLS_SAVE_THREAD_INFO(prefix) do { \
    SF_scal_save(prefix "_threads_max", (double)omp_get_max_threads()); \
    SF_scal_save(prefix "_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED); \
} while(0)

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
   I/O and CSV Processing Constants
   ============================================================================ */

/*
    I/O buffer size for file operations.
    Larger buffers reduce syscalls but use more memory.
    64KB is optimal for most systems.
*/
#define CTOOLS_IO_BUFFER_SIZE    (64 * 1024)

/*
    Chunk sizes for parallel CSV processing.
    CTOOLS_IMPORT_CHUNK_SIZE: bytes per parsing chunk (8MB default)
    CTOOLS_EXPORT_CHUNK_SIZE: rows per formatting chunk (10K default)
*/
#define CTOOLS_IMPORT_CHUNK_SIZE (8 * 1024 * 1024)
#define CTOOLS_EXPORT_CHUNK_SIZE 10000

/*
    Maximum threads for I/O operations.
    May differ from NUM_THREADS for I/O-bound vs CPU-bound work.
*/
#define CTOOLS_IO_MAX_THREADS    12

/*
    Arena allocator block size for string pooling.
    1MB blocks reduce allocation overhead for many small strings.
*/
#define CTOOLS_ARENA_BLOCK_SIZE  (1024 * 1024)

/*
    Stata variable limits.
*/
#define CTOOLS_MAX_VARNAME_LEN   32
#define CTOOLS_MAX_STRING_LEN    2045
#define CTOOLS_MAX_COLUMNS       32767

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

/* ============================================================================
   Safe size arithmetic for memory allocation

   Prevents integer overflow when computing allocation sizes like N*K*sizeof(T).
   Returns 0 on success, -1 on overflow.
   ============================================================================ */

#include <stdint.h>

/*
    Multiply two size_t values safely, checking for overflow.
    Returns 0 on success (result stored in *out), -1 on overflow.
*/
static inline int ctools_safe_mul_size(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) {
        *out = 0;
        return -1;  /* overflow */
    }
    *out = a * b;
    return 0;
}

/*
    Compute a * b * c safely for allocation sizes.
    Typical use: ctools_safe_alloc_size(N, K, sizeof(double), &size)
    Returns 0 on success, -1 on overflow.
*/
static inline int ctools_safe_alloc_size(size_t a, size_t b, size_t c, size_t *out)
{
    size_t tmp;
    if (ctools_safe_mul_size(a, b, &tmp) != 0) {
        *out = 0;
        return -1;
    }
    if (ctools_safe_mul_size(tmp, c, out) != 0) {
        *out = 0;
        return -1;
    }
    return 0;
}

/*
    Allocate memory with overflow-checked size computation.
    Returns NULL on overflow or allocation failure.
    Usage: ptr = ctools_safe_malloc2(N, sizeof(double));
           ptr = ctools_safe_malloc3(N, K, sizeof(double));
*/
static inline void *ctools_safe_malloc2(size_t a, size_t b)
{
    size_t size;
    if (ctools_safe_mul_size(a, b, &size) != 0) return NULL;
    return malloc(size);
}

static inline void *ctools_safe_malloc3(size_t a, size_t b, size_t c)
{
    size_t size;
    if (ctools_safe_alloc_size(a, b, c, &size) != 0) return NULL;
    return malloc(size);
}

static inline void *ctools_safe_calloc2(size_t a, size_t b)
{
    size_t size;
    if (ctools_safe_mul_size(a, b, &size) != 0) return NULL;
    return calloc(1, size);
}

static inline void *ctools_safe_calloc3(size_t a, size_t b, size_t c)
{
    size_t size;
    if (ctools_safe_alloc_size(a, b, c, &size) != 0) return NULL;
    return calloc(1, size);
}

/*
    Overflow-safe aligned allocation for (count * element_size).
    Returns NULL on overflow or allocation failure.
    MUST be freed with ctools_aligned_free().
*/
static inline void *ctools_safe_aligned_alloc2(size_t alignment, size_t count, size_t element_size)
{
    size_t size;
    if (ctools_safe_mul_size(count, element_size, &size) != 0) {
        return NULL;  /* Overflow */
    }
    return ctools_aligned_alloc(alignment, size);
}

#endif /* CTOOLS_CONFIG_H */

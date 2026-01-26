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
    Minimum observations per thread for parallel sort.
    Below this, sequential radix sort is faster due to better cache locality.
*/
#define MIN_OBS_PER_THREAD 100000

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

/* ============================================================================
   Adaptive Prefetch Distances

   Prefetch distances are scaled based on detected cache hierarchy.
   Larger L1 caches allow larger prefetch distances without evicting
   useful data. These values are computed once at startup.
   ============================================================================ */

/* Platform-specific cache detection headers */
#if defined(__APPLE__)
    #include <sys/sysctl.h>
#elif defined(__linux__)
    #include <unistd.h>
#endif

/*
    Cache size detection structure.
    Sizes are in bytes, 0 means unknown/undetected.
*/
typedef struct {
    size_t l1d_size;    /* L1 data cache size */
    size_t l2_size;     /* L2 cache size */
    size_t l3_size;     /* L3 cache size */
    size_t line_size;   /* Cache line size */
} ctools_cache_info;

/*
    Detect cache sizes at runtime.
    Returns cache info structure with detected sizes.
    Unknown values are set to reasonable defaults.
*/
static inline ctools_cache_info ctools_detect_cache_sizes(void)
{
    ctools_cache_info info = {0, 0, 0, CACHE_LINE_SIZE};

#if defined(__APPLE__)
    /* macOS: use sysctl */
    size_t size = sizeof(size_t);

    /* L1 data cache */
    if (sysctlbyname("hw.l1dcachesize", &info.l1d_size, &size, NULL, 0) != 0) {
        info.l1d_size = 32 * 1024;  /* Default 32KB */
    }

    /* L2 cache */
    size = sizeof(size_t);
    if (sysctlbyname("hw.l2cachesize", &info.l2_size, &size, NULL, 0) != 0) {
        info.l2_size = 256 * 1024;  /* Default 256KB */
    }

    /* L3 cache (may not exist on all systems) */
    size = sizeof(size_t);
    if (sysctlbyname("hw.l3cachesize", &info.l3_size, &size, NULL, 0) != 0) {
        info.l3_size = 0;  /* Not available */
    }

    /* Cache line size */
    size = sizeof(size_t);
    if (sysctlbyname("hw.cachelinesize", &info.line_size, &size, NULL, 0) != 0) {
        info.line_size = 64;
    }

#elif defined(__linux__)
    /* Linux: use sysconf */
    long l1d = sysconf(_SC_LEVEL1_DCACHE_SIZE);
    long l2 = sysconf(_SC_LEVEL2_CACHE_SIZE);
    long l3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
    long line = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

    info.l1d_size = (l1d > 0) ? (size_t)l1d : 32 * 1024;
    info.l2_size = (l2 > 0) ? (size_t)l2 : 256 * 1024;
    info.l3_size = (l3 > 0) ? (size_t)l3 : 0;
    info.line_size = (line > 0) ? (size_t)line : 64;

#else
    /* Fallback: use conservative defaults */
    info.l1d_size = 32 * 1024;   /* 32KB L1 */
    info.l2_size = 256 * 1024;   /* 256KB L2 */
    info.l3_size = 0;            /* Unknown L3 */
    info.line_size = 64;
#endif

    return info;
}

/*
    Adaptive prefetch distances based on cache sizes.

    The idea: larger L1 cache allows prefetching more data ahead
    without evicting useful data. We scale prefetch distances based
    on detected cache sizes.

    Base distances (for 32KB L1):
    - NEAR:   8 elements  (64 bytes = 1 cache line)
    - FAR:   32 elements (256 bytes = 4 cache lines)
    - STREAM: 64 elements (512 bytes = 8 cache lines)

    For larger caches, we can increase these proportionally.
*/
typedef struct {
    int near_dist;  /* Near prefetch distance (into L1) */
    int far_dist;   /* Far prefetch distance (into L2/L3) */
    int stream;     /* Streaming prefetch distance */
    int general;    /* General prefetch distance */
} ctools_prefetch_distances;

/*
    Compute optimal prefetch distances based on cache sizes.
    Call this once at startup and cache the result.
*/
static inline ctools_prefetch_distances ctools_compute_prefetch_distances(void)
{
    ctools_cache_info cache = ctools_detect_cache_sizes();
    ctools_prefetch_distances dist;

    /* Scale factor based on L1 cache size (baseline: 32KB) */
    /* Larger L1 = can prefetch more ahead without eviction */
    int scale = 1;
    if (cache.l1d_size >= 128 * 1024) {
        scale = 4;  /* 128KB+ L1 (Apple M1/M2 have 192KB) */
    } else if (cache.l1d_size >= 64 * 1024) {
        scale = 2;  /* 64KB L1 */
    }

    /* Base distances scaled by cache size */
    dist.near_dist = 8 * scale;      /* L1 prefetch: 8-32 elements */
    dist.far_dist = 32 * scale;      /* L2/L3 prefetch: 32-128 elements */
    dist.stream = 64 * scale;        /* Streaming: 64-256 elements */
    dist.general = 16 * scale;       /* General purpose: 16-64 elements */

    /* Cap at reasonable maximums to avoid prefetching too far ahead */
    if (dist.near_dist > 64) dist.near_dist = 64;
    if (dist.far_dist > 256) dist.far_dist = 256;
    if (dist.stream > 512) dist.stream = 512;
    if (dist.general > 128) dist.general = 128;

    return dist;
}

/*
    Global prefetch distances - initialized once.
    Use ctools_get_prefetch_distances() to access.

    Thread-safety: Uses atomic operations for the initialization flag to prevent
    data races. Multiple threads may compute distances simultaneously during
    first access (benign race - all produce same result), but the flag ensures
    we don't read partially-written struct fields.
*/
static ctools_prefetch_distances _ctools_prefetch_dist = {0, 0, 0, 0};
#if defined(__GNUC__) || defined(__clang__)
static volatile int _ctools_prefetch_initialized = 0;
#elif defined(_WIN32) && defined(_MSC_VER)
static volatile long _ctools_prefetch_initialized = 0;  /* LONG for Interlocked* */
#else
static int _ctools_prefetch_initialized = 0;
#endif

/*
    Get adaptive prefetch distances (lazy initialization).
    Thread-safe using double-checked locking pattern with memory barriers.
*/
static inline ctools_prefetch_distances ctools_get_prefetch_distances(void)
{
    /* Fast path with acquire semantics */
#if defined(__GNUC__) || defined(__clang__)
    if (__atomic_load_n(&_ctools_prefetch_initialized, __ATOMIC_ACQUIRE)) {
        return _ctools_prefetch_dist;
    }
#elif defined(_WIN32) && defined(_MSC_VER)
    /* Use InterlockedOr for atomic read with acquire semantics.
     * Prevents torn reads on ARM64 Windows. */
    if (InterlockedOr(&_ctools_prefetch_initialized, 0)) {
        return _ctools_prefetch_dist;
    }
#else
    if (_ctools_prefetch_initialized) {
        return _ctools_prefetch_dist;
    }
#endif

    /* Slow path: compute distances.
     * Multiple threads may compute simultaneously - this is a benign race
     * since all threads produce identical results. The atomic store ensures
     * subsequent readers see fully-written struct. */
    ctools_prefetch_distances dist = ctools_compute_prefetch_distances();
    _ctools_prefetch_dist = dist;

    /* Release semantics: ensure struct write is visible before flag */
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&_ctools_prefetch_initialized, 1, __ATOMIC_RELEASE);
#elif defined(_WIN32) && defined(_MSC_VER)
    /* InterlockedExchange provides release semantics on all Windows architectures */
    InterlockedExchange(&_ctools_prefetch_initialized, 1);
#else
    _ctools_prefetch_initialized = 1;
#endif

    return dist;
}

/*
    Legacy compatibility: fixed prefetch distance macro.
    New code should use ctools_get_prefetch_distances() instead.
*/
#define PREFETCH_DISTANCE 16

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
    SF_scal_save(prefix "_threads_max", (double)ctools_get_max_threads()); \
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

/* ============================================================================
   Unified Thread Management

   All ctools modules use these functions for thread configuration.
   The thread limit is determined at runtime using omp_get_max_threads(),
   and can be overridden via the threads() option in Stata commands.

   Thread limit hierarchy:
   1. User-specified via threads() option -> ctools_set_max_threads()
   2. omp_get_max_threads() (respects OMP_NUM_THREADS env var)
   3. Falls back to 1 if OpenMP unavailable

   Usage:
   - Call ctools_get_max_threads() to get current limit
   - Call ctools_set_max_threads(n) to override (0 = reset to default)
   ============================================================================ */

/* Get the current maximum thread count for ctools operations.
 * Returns the user-set limit if specified, otherwise omp_get_max_threads(). */
int ctools_get_max_threads(void);

/* Set the maximum thread count for ctools operations.
 * Pass 0 to reset to default (omp_get_max_threads()).
 * Pass n > 0 to set a specific limit. */
void ctools_set_max_threads(int n);

/* Reset thread limit to default (runtime-detected) */
void ctools_reset_max_threads(void);

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
    #define ctools_memory_barrier() do { _ReadWriteBarrier(); MemoryBarrier(); } while(0)
#elif defined(__GNUC__) || defined(__clang__)
    #define ctools_memory_barrier() __sync_synchronize()
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_ATOMICS__)
    /* C11 atomics available */
    #include <stdatomic.h>
    #define ctools_memory_barrier() atomic_thread_fence(memory_order_seq_cst)
#else
    /* Fallback: volatile access acts as a compiler barrier.
     * This provides no hardware memory ordering guarantees but prevents
     * compiler reordering. Safe for single-threaded or already-locked contexts. */
    static volatile int _ctools_barrier_dummy = 0;
    #define ctools_memory_barrier() do { _ctools_barrier_dummy = _ctools_barrier_dummy; } while(0)
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

/* ============================================================================
   Convenience macros for cache-line aligned allocation

   These macros provide a concise way to allocate memory aligned to the
   cache line boundary, which is the most common alignment requirement.
   ============================================================================ */

/*
    Allocate memory aligned to cache line boundary.
    Returns NULL on failure.
    MUST be freed with ctools_aligned_free().
*/
#define ctools_cacheline_alloc(size) \
    ctools_aligned_alloc(CACHE_LINE_SIZE, (size))

/*
    Overflow-safe cache-line aligned allocation for (count * element_size).
    Returns NULL on overflow or allocation failure.
    MUST be freed with ctools_aligned_free().
*/
#define ctools_safe_cacheline_alloc2(count, elem_size) \
    ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, (count), (elem_size))

#endif /* CTOOLS_CONFIG_H */

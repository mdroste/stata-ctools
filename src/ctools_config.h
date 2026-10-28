/*
    ctools_config.h
    Configuration parameters for ctools performance tuning

    These values can be adjusted based on profiling results and target hardware.
*/

#ifndef CTOOLS_CONFIG_H
#define CTOOLS_CONFIG_H

#include <stddef.h>  /* size_t */

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

#include "ctools_platform.h"

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

#include "ctools_alloc.h"

/* ============================================================================
   I/O and CSV Processing Constants
   ============================================================================ */

/*
    I/O buffer size for file operations.
    Larger buffers reduce syscalls but use more memory.
    64KB is optimal for most systems.
*/
#define CTOOLS_IO_BUFFER_SIZE    (64 * 1024)

/* Rows per formatting chunk for parallel CSV export. */
#define CTOOLS_EXPORT_CHUNK_SIZE 10000

/* ============================================================================
   Unified Thread Management

   All ctools modules use these functions for thread configuration.
   The default thread limit is resolved once per process at plugin entry,
   and can be overridden via the threads() option in Stata commands.

   Thread limit hierarchy (first that is set wins):
   1. threads() option on the command            -> ctools_set_max_threads()
   2. Session setting: global CTOOLS_THREADS     (ctools, threads(#))
   3. CTOOLS_NUM_THREADS environment variable
   4. Saved per-user setting                     (ctools, threads(#) permanently)
   5. OMP_NUM_THREADS environment variable (first entry)
   6. CTOOLS_DEFAULT_THREAD_CAP
   Every source is clamped to the CPUs this process may use (affinity mask and
   Linux cgroup quota), OMP_THREAD_LIMIT and CTOOLS_MAX_THREADS. The default
   is captured before third-party plugins can alter the OpenMP runtime.
   OpenMP-only kernels use one worker when OpenMP is unavailable.

   Usage:
   - Call ctools_get_max_threads() to get current limit
   - Call ctools_set_max_threads(n) to override (0 = reset to default)
   ============================================================================ */

/* Default thread count when nothing is configured. Most ctools work is memory
 * bound, so large servers gain little beyond this and are often shared. */
#define CTOOLS_DEFAULT_THREAD_CAP 16

typedef struct {
    int detected;    /* logical CPUs reported by the OS */
    int available;   /* CPUs this process may use: affinity, cgroup quota */
    int omp_limit;   /* OMP_THREAD_LIMIT: >0 set, 0 unset, -1 invalid */
    int env_ctools;  /* CTOOLS_NUM_THREADS, same encoding */
    int saved;       /* threads= in the settings file, same encoding */
    int env_omp;     /* OMP_NUM_THREADS (first entry), same encoding */
    int limit;       /* most threads any source may request */
    int threads;     /* resolved default */
    const char *source;  /* which input chose threads */
    char settings_path[4096];  /* UTF-8; empty when no home directory */
} ctools_thread_policy;

/* Parse a thread count: >0 valid, 0 empty, -1 invalid. When list is nonzero
 * the value may continue with ",..." (OMP_NUM_THREADS nesting levels). */
int ctools_parse_thread_value(const char *text, int list);

/* Fill limit, threads and source from the other fields. Pure. */
void ctools_thread_policy_resolve(ctools_thread_policy *policy);

/* Process-wide policy, read from the OS and environment on first use. */
const ctools_thread_policy *ctools_thread_policy_get(void);

/* Clamp an explicit request to the policy limit. */
int ctools_clamp_thread_request(int n);

/* Save (n > 0) or remove (n == 0) the per-user thread setting and make it
 * the default for this process. Returns 0, or -1 with errno set. */
int ctools_save_thread_setting(int n);

/* Capture the resolved default at plugin init, before external
 * plugins can pollute it. Call once at plugin entry. */
void ctools_init_thread_count(void);

/* Get the current maximum thread count for ctools operations.
 * Returns the user-set limit if specified, otherwise the default
 * resolved at init (immune to external OpenMP pollution). */
/* Bound explicit requests before passing them to OpenMP. Large scratch users
 * further reduce concurrency to a 256 MiB aggregate workspace target. One
 * worker is always allowed, even when its required buffer exceeds the target. */
#define CTOOLS_MAX_THREADS 256
#define CTOOLS_WORKSPACE_BUDGET ((size_t)256 * 1024 * 1024)
static inline int ctools_workspace_threads(int requested, size_t tasks, size_t bytes_per_worker)
{
    size_t n = requested > 0 ? (size_t)requested : 1;
    if (n > CTOOLS_MAX_THREADS) n = CTOOLS_MAX_THREADS;
    if (tasks && n > tasks) n = tasks;
    if (bytes_per_worker) {
        size_t affordable = CTOOLS_WORKSPACE_BUDGET / bytes_per_worker;
        if (!affordable) affordable = 1;
        if (n > affordable) n = affordable;
    }
    return (int)n;
}

int ctools_get_max_threads(void);

/* OpenMP kernels must not confuse pthread/CPU capacity with an OMP team.
 * Manual partitions still need work-sharing loops: a runtime may supply fewer
 * workers than requested, including inside nested parallel regions. */
static inline int ctools_get_openmp_threads(void)
{
#ifdef _OPENMP
    return ctools_get_max_threads();
#else
    return 1;
#endif
}

/* Set the maximum thread count for ctools operations.
 * Pass 0 to reset to the captured system default.
 * Pass n > 0 to set a specific limit. */
void ctools_set_max_threads(int n);

/* Reset thread limit to the resolved default.
 * Also restores the OpenMP runtime to that default. */
void ctools_reset_max_threads(void);

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

#endif /* CTOOLS_CONFIG_H */

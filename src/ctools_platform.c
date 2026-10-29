#include "ctools_config.h"
#include <stdatomic.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

uint64_t ctools_physical_memory_bytes(void)
{
#if defined(__APPLE__)
    uint64_t bytes = 0;
    size_t size = sizeof(bytes);
    return sysctlbyname("hw.memsize", &bytes, &size, NULL, 0) == 0 ? bytes : 0;
#elif defined(_WIN32)
    MEMORYSTATUSEX status = {0};
    status.dwLength = sizeof(status);
    return GlobalMemoryStatusEx(&status) ? status.ullTotalPhys : 0;
#elif defined(__linux__)
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0 && (uint64_t)pages <= UINT64_MAX / (uint64_t)page_size)
        return (uint64_t)pages * (uint64_t)page_size;
#endif
    return 0;
}

ctools_cache_info ctools_detect_cache_sizes(void)
{
    ctools_cache_info info = {0, 0, 0, CACHE_LINE_SIZE};

#if defined(__APPLE__)
    /* macOS: use sysctl */
    size_t size = sizeof(size_t);

    /* L1 data cache */
    if (sysctlbyname("hw.l1dcachesize", &info.l1d_size, &size, NULL, 0) != 0)
    {
        info.l1d_size = 32 * 1024; /* Default 32KB */
    }

    /* L2 cache */
    size = sizeof(size_t);
    if (sysctlbyname("hw.l2cachesize", &info.l2_size, &size, NULL, 0) != 0)
    {
        info.l2_size = 256 * 1024; /* Default 256KB */
    }

    /* L3 cache (may not exist on all systems) */
    size = sizeof(size_t);
    if (sysctlbyname("hw.l3cachesize", &info.l3_size, &size, NULL, 0) != 0)
    {
        info.l3_size = 0; /* Not available */
    }

    /* Cache line size */
    size = sizeof(size_t);
    if (sysctlbyname("hw.cachelinesize", &info.line_size, &size, NULL, 0) != 0)
    {
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
    info.l1d_size = 32 * 1024; /* 32KB L1 */
    info.l2_size = 256 * 1024; /* 256KB L2 */
    info.l3_size = 0;          /* Unknown L3 */
    info.line_size = 64;
#endif

    return info;
}

ctools_prefetch_distances ctools_compute_prefetch_distances(void)
{
    ctools_cache_info cache = ctools_detect_cache_sizes();
    ctools_prefetch_distances dist;

    /* Scale factor based on L1 cache size (baseline: 32KB) */
    /* Larger L1 = can prefetch more ahead without eviction */
    int scale = 1;
    if (cache.l1d_size >= 128 * 1024)
    {
        scale = 4; /* 128KB+ L1 (Apple M1/M2 have 192KB) */
    }
    else if (cache.l1d_size >= 64 * 1024)
    {
        scale = 2; /* 64KB L1 */
    }

    /* Base distances scaled by cache size */
    dist.near_dist = 8 * scale; /* L1 prefetch: 8-32 elements */
    dist.far_dist = 32 * scale; /* L2/L3 prefetch: 32-128 elements */
    dist.stream = 64 * scale;   /* Streaming: 64-256 elements */
    dist.general = 16 * scale;  /* General purpose: 16-64 elements */

    /* Cap at reasonable maximums to avoid prefetching too far ahead */
    if (dist.near_dist > 64)
        dist.near_dist = 64;
    if (dist.far_dist > 256)
        dist.far_dist = 256;
    if (dist.stream > 512)
        dist.stream = 512;
    if (dist.general > 128)
        dist.general = 128;

    return dist;
}

static atomic_int prefetch_state = 0;
static ctools_prefetch_distances prefetch_distances;

ctools_prefetch_distances ctools_get_prefetch_distances(void)
{
    if (atomic_load_explicit(&prefetch_state, memory_order_acquire) == 2)
        return prefetch_distances;
    int expected = 0;
    if (atomic_compare_exchange_strong_explicit(&prefetch_state, &expected, 1, memory_order_acq_rel,
                                                memory_order_acquire))
    {
        prefetch_distances = ctools_compute_prefetch_distances();
        atomic_store_explicit(&prefetch_state, 2, memory_order_release);
    }
    else
    {
        while (atomic_load_explicit(&prefetch_state, memory_order_acquire) != 2)
        {
        }
    }
    return prefetch_distances;
}

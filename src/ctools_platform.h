#ifndef CTOOLS_PLATFORM_H
#define CTOOLS_PLATFORM_H
#include <stddef.h>
typedef struct
{
    size_t l1d_size, l2_size, l3_size, line_size;
} ctools_cache_info;
typedef struct
{
    int near_dist, far_dist, stream, general;
} ctools_prefetch_distances;
ctools_cache_info ctools_detect_cache_sizes(void);
ctools_prefetch_distances ctools_compute_prefetch_distances(void);
/* One process-wide immutable snapshot, published by a single initializer. */
ctools_prefetch_distances ctools_get_prefetch_distances(void);
#endif

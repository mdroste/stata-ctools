#ifndef CTOOLS_SORT_BLOCK_H
#define CTOOLS_SORT_BLOCK_H
#include "ctools_data.h"
#include "ctools_config.h"
#include <string.h>
/* Stable eight-pass LSD block primitive. Keys indexed by order; scratch holds
 * len indices and must not alias order. Thresholds belong to each caller. */
static inline void ctools_sort_block_radix(perm_idx_t *CTOOLS_RESTRICT order,
                                           const uint64_t *CTOOLS_RESTRICT keys, size_t start,
                                           size_t len, perm_idx_t *CTOOLS_RESTRICT temp)
{
    size_t counts[256];
    size_t offsets[256];
    perm_idx_t *src = order + start;
    perm_idx_t *dst = temp;

    for (int byte_pos = 0; byte_pos < 8; byte_pos++)
    {
        int shift = byte_pos * 8;

        memset(counts, 0, sizeof(counts));

        /* Counting with prefetch and unrolling */
        size_t i = 0;
        for (; i + 4 <= len; i += 4)
        {
            if (i + PREFETCH_DISTANCE < len)
            {
                CTOOLS_PREFETCH(&keys[src[i + PREFETCH_DISTANCE]]);
            }
            counts[(keys[src[i]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 1]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 2]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 3]] >> shift) & 0xFF]++;
        }
        for (; i < len; i++)
        {
            counts[(keys[src[i]] >> shift) & 0xFF]++;
        }

        /* Check for uniform distribution */
        int non_empty = 0;
        for (size_t b = 0; b < 256; b++)
        {
            if (counts[b] > 0)
            {
                non_empty++;
                if (non_empty > 1)
                    break;
            }
        }
        if (non_empty <= 1)
            continue;

        /* Prefix sum */
        offsets[0] = 0;
        for (size_t b = 1; b < 256; b++)
        {
            offsets[b] = offsets[b - 1] + counts[b - 1];
        }

        /* Scatter with prefetch */
        for (i = 0; i < len; i++)
        {
            if (i + PREFETCH_DISTANCE < len)
            {
                CTOOLS_PREFETCH(&keys[src[i + PREFETCH_DISTANCE]]);
            }
            uint8_t byte_val = (keys[src[i]] >> shift) & 0xFF;
            dst[offsets[byte_val]++] = src[i];
        }

        /* Swap */
        perm_idx_t *tmp = src;
        src = dst;
        dst = tmp;
    }

    /* Copy back if needed */
    if (src != order + start)
    {
        memcpy(order + start, src, len * sizeof(perm_idx_t));
    }
}
#endif

/*
 * cmerge_radix_sort.c
 * Parallel LSD Radix Sort for order pairs
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "cmerge_radix_sort.h"
#include "../ctools_threads.h"
#include "../ctools_config.h"

/* Thread arguments for parallel histogram */
typedef struct {
    cmerge_order_pair_t *pairs;
    size_t start;
    size_t end;
    int shift;
    size_t *local_counts;  /* [256] */
} cmerge_hist_args_t;

/* Thread arguments for parallel scatter */
typedef struct {
    cmerge_order_pair_t *src;
    cmerge_order_pair_t *dst;
    size_t start;
    size_t end;
    int shift;
    size_t *offsets;  /* Starting offset for each bucket for this thread */
} cmerge_scatter_args_t;

/* Comparison function for qsort - makes sort stable by comparing orig_idx on ties */
int cmerge_compare_order_pairs(const void *a, const void *b)
{
    const cmerge_order_pair_t *pa = (const cmerge_order_pair_t *)a;
    const cmerge_order_pair_t *pb = (const cmerge_order_pair_t *)b;
    if (pa->order_key < pb->order_key) return -1;
    if (pa->order_key > pb->order_key) return 1;
    /* Stable sort: preserve original order for equal keys */
    if (pa->orig_idx < pb->orig_idx) return -1;
    if (pa->orig_idx > pb->orig_idx) return 1;
    return 0;
}

static void *cmerge_histogram_thread(void *arg)
{
    cmerge_hist_args_t *args = (cmerge_hist_args_t *)arg;
    cmerge_order_pair_t *pairs = args->pairs;
    int shift = args->shift;
    size_t *counts = args->local_counts;

    memset(counts, 0, 256 * sizeof(size_t));

    for (size_t i = args->start; i < args->end; i++) {
        uint8_t digit = (pairs[i].order_key >> shift) & 0xFF;
        counts[digit]++;
    }

    return NULL;
}

static void *cmerge_scatter_thread(void *arg)
{
    cmerge_scatter_args_t *args = (cmerge_scatter_args_t *)arg;
    cmerge_order_pair_t *src = args->src;
    cmerge_order_pair_t *dst = args->dst;
    int shift = args->shift;
    size_t *offsets = args->offsets;

    for (size_t i = args->start; i < args->end; i++) {
        uint8_t digit = (src[i].order_key >> shift) & 0xFF;
        dst[offsets[digit]++] = src[i];
    }

    return NULL;
}

void cmerge_radix_sort_order_pairs(cmerge_order_pair_t *pairs, size_t n)
{
    /* For small arrays, qsort is faster due to parallel overhead */
    if (n < 10000) {
        qsort(pairs, n, sizeof(cmerge_order_pair_t), cmerge_compare_order_pairs);
        return;
    }

    /* Determine thread count */
    int num_threads = ctools_get_max_threads();
    if ((size_t)num_threads > n / 1000) num_threads = (int)(n / 1000);
    if (num_threads < 1) num_threads = 1;

    /* Get persistent thread pool */
    ctools_persistent_pool *pool = ctools_get_global_pool();

    /* Allocate auxiliary buffer and thread resources (overflow-safe) */
    cmerge_order_pair_t *aux = (cmerge_order_pair_t *)ctools_safe_malloc2(n, sizeof(cmerge_order_pair_t));
    cmerge_hist_args_t *hist_args = (cmerge_hist_args_t *)ctools_safe_malloc2(num_threads, sizeof(cmerge_hist_args_t));
    cmerge_scatter_args_t *scatter_args = (cmerge_scatter_args_t *)ctools_safe_malloc2(num_threads, sizeof(cmerge_scatter_args_t));
    size_t *all_counts = (size_t *)ctools_safe_malloc3(num_threads, 256, sizeof(size_t));
    size_t *all_offsets = (size_t *)ctools_safe_malloc3(num_threads, 256, sizeof(size_t));

    if (!aux || !hist_args || !scatter_args || !all_counts || !all_offsets) {
        /* Fallback to qsort on allocation failure */
        free(aux);
        free(hist_args);
        free(scatter_args);
        free(all_counts);
        free(all_offsets);
        qsort(pairs, n, sizeof(cmerge_order_pair_t), cmerge_compare_order_pairs);
        return;
    }

    /* Calculate chunk boundaries */
    size_t chunk_size = n / num_threads;
    size_t remainder = n % num_threads;

    cmerge_order_pair_t *src = pairs;
    cmerge_order_pair_t *dst = aux;

    /* 4 passes of 8-bit radix sort on order_key (32 bits sufficient for Stata) */
    for (int pass = 0; pass < 4; pass++) {
        int shift = pass * 8;

        /* Phase 1: Parallel histogram */
        size_t offset = 0;
        for (int t = 0; t < num_threads; t++) {
            hist_args[t].pairs = src;
            hist_args[t].start = offset;
            hist_args[t].end = offset + chunk_size + (t < (int)remainder ? 1 : 0);
            hist_args[t].shift = shift;
            hist_args[t].local_counts = all_counts + t * 256;
            offset = hist_args[t].end;
        }

        if (pool != NULL && num_threads >= 2) {
            ctools_persistent_pool_submit_batch(pool, cmerge_histogram_thread,
                                                 hist_args, num_threads,
                                                 sizeof(cmerge_hist_args_t));
            ctools_persistent_pool_wait(pool);
        } else {
            for (int t = 0; t < num_threads; t++) {
                cmerge_histogram_thread(&hist_args[t]);
            }
        }

        /* Compute global counts and per-thread offsets */
        size_t global_counts[256];
        memset(global_counts, 0, sizeof(global_counts));

        for (int t = 0; t < num_threads; t++) {
            size_t *local = all_counts + t * 256;
            for (int b = 0; b < 256; b++) {
                global_counts[b] += local[b];
            }
        }

        /* Convert global counts to prefix sums */
        size_t total = 0;
        size_t global_offsets[256];
        for (int b = 0; b < 256; b++) {
            global_offsets[b] = total;
            total += global_counts[b];
        }

        /* Compute per-thread starting offsets for each bucket */
        for (int b = 0; b < 256; b++) {
            size_t bucket_offset = global_offsets[b];
            for (int t = 0; t < num_threads; t++) {
                all_offsets[t * 256 + b] = bucket_offset;
                bucket_offset += all_counts[t * 256 + b];
            }
        }

        /* Phase 2: Parallel scatter */
        offset = 0;
        for (int t = 0; t < num_threads; t++) {
            scatter_args[t].src = src;
            scatter_args[t].dst = dst;
            scatter_args[t].start = offset;
            scatter_args[t].end = offset + chunk_size + (t < (int)remainder ? 1 : 0);
            scatter_args[t].shift = shift;
            scatter_args[t].offsets = all_offsets + t * 256;
            offset = scatter_args[t].end;
        }

        if (pool != NULL && num_threads >= 2) {
            ctools_persistent_pool_submit_batch(pool, cmerge_scatter_thread,
                                                 scatter_args, num_threads,
                                                 sizeof(cmerge_scatter_args_t));
            ctools_persistent_pool_wait(pool);
        } else {
            for (int t = 0; t < num_threads; t++) {
                cmerge_scatter_thread(&scatter_args[t]);
            }
        }

        /* Swap buffers */
        cmerge_order_pair_t *tmp = src;
        src = dst;
        dst = tmp;
    }

    /* If result is in aux, copy back to pairs */
    if (src != pairs) {
        memcpy(pairs, src, n * sizeof(cmerge_order_pair_t));
    }

    free(aux);
    free(hist_args);
    free(scatter_args);
    free(all_counts);
    free(all_offsets);
}

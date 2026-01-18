/*
 * cbinscatter_bins.c
 *
 * Bin computation implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 *
 * Optimizations:
 * - Direct bin assignment from sorted order (no binary search)
 * - O(N) weighted quantile computation via cumulative weights
 * - Single-pass bin statistics accumulation
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "stplugin.h"
#include "cbinscatter_bins.h"
#include "../ctools_config.h"

/* Stata missing value check */
#define STATA_MISSING 8.988465674311579e+307
#define IS_MISSING(x) ((x) >= STATA_MISSING)

/* Threshold for using histogram-based binning vs full sort */
#define HISTOGRAM_THRESHOLD 50000
#define HISTOGRAM_BUCKETS 4096

/* ========================================================================
 * Quickselect for finding k-th element
 * ======================================================================== */

static inline __attribute__((unused)) void swap_double(ST_double *a, ST_double *b) {
    ST_double t = *a; *a = *b; *b = t;
}

static inline void swap_int(ST_int *a, ST_int *b) {
    ST_int t = *a; *a = *b; *b = t;
}

/*
 * Partition array around pivot, return pivot position
 * Uses median-of-three pivot selection for better performance
 */
static ST_int partition_indexed(ST_double *vals, ST_int *idx, ST_int lo, ST_int hi) {
    /* Median-of-three pivot selection */
    ST_int mid = lo + (hi - lo) / 2;

    /* Sort lo, mid, hi */
    if (vals[idx[mid]] < vals[idx[lo]]) { swap_int(&idx[lo], &idx[mid]); }
    if (vals[idx[hi]] < vals[idx[lo]]) { swap_int(&idx[lo], &idx[hi]); }
    if (vals[idx[hi]] < vals[idx[mid]]) { swap_int(&idx[mid], &idx[hi]); }

    /* Use median as pivot, move to hi-1 */
    swap_int(&idx[mid], &idx[hi - 1]);
    ST_double pivot = vals[idx[hi - 1]];

    ST_int i = lo;
    ST_int j = hi - 1;

    for (;;) {
        while (vals[idx[++i]] < pivot) {}
        while (vals[idx[--j]] > pivot) { if (j == lo) break; }
        if (i >= j) break;
        swap_int(&idx[i], &idx[j]);
    }

    swap_int(&idx[i], &idx[hi - 1]);
    return i;
}

/*
 * Quickselect: find k-th smallest element
 * After this call, idx[k] contains the index of the k-th smallest value
 * and all elements before k are <= vals[idx[k]]
 */
static __attribute__((unused)) void quickselect_indexed(ST_double *vals, ST_int *idx, ST_int lo, ST_int hi, ST_int k) {
    while (hi > lo) {
        if (hi - lo < 10) {
            /* Insertion sort for small arrays */
            for (ST_int i = lo + 1; i <= hi; i++) {
                ST_int temp = idx[i];
                ST_double val = vals[temp];
                ST_int j = i;
                while (j > lo && vals[idx[j - 1]] > val) {
                    idx[j] = idx[j - 1];
                    j--;
                }
                idx[j] = temp;
            }
            return;
        }

        ST_int pivot_pos = partition_indexed(vals, idx, lo, hi);

        if (k == pivot_pos) {
            return;
        } else if (k < pivot_pos) {
            hi = pivot_pos - 1;
        } else {
            lo = pivot_pos + 1;
        }
    }
}

/* ========================================================================
 * Argsort comparison (for fallback full sort)
 * ======================================================================== */

/* Global pointer for qsort comparison (thread-local would be better but this is simpler) */
static const ST_double *g_sort_values;

static int compare_indices(const void *a, const void *b) {
    ST_int ia = *(const ST_int *)a;
    ST_int ib = *(const ST_int *)b;
    ST_double va = g_sort_values[ia];
    ST_double vb = g_sort_values[ib];

    /* Handle missing values - sort to end */
    if (IS_MISSING(va) && IS_MISSING(vb)) return 0;
    if (IS_MISSING(va)) return 1;
    if (IS_MISSING(vb)) return -1;

    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

/* ========================================================================
 * Bin Statistics (optimized single-pass)
 * ======================================================================== */

ST_retcode compute_bin_statistics(
    const ST_double *y,
    const ST_double *x,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int nquantiles,
    ST_int compute_se,
    ByGroupResult *result
) {
    ST_int i, b;
    ST_double *sum_y = NULL, *sum_x = NULL;
    ST_double *sum_y2 = NULL, *sum_x2 = NULL;
    ST_double *sum_w = NULL;
    ST_int *counts = NULL;
    ST_double w, yi, xi, mean_y, mean_x, var_y, var_x;
    ST_int actual_bins = 0;

    /* Validate nquantiles to prevent overflow */
    if (nquantiles <= 0) {
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Allocate accumulators - all in one block for cache efficiency */
    /* Use safe multiplication to prevent overflow */
    size_t base_per_bin = 3 * sizeof(ST_double) + sizeof(ST_int);
    size_t alloc_size;
    if (ctools_safe_mul_size((size_t)nquantiles, base_per_bin, &alloc_size) != 0) {
        return CBINSCATTER_ERR_MEMORY;
    }
    if (compute_se) {
        size_t se_size;
        if (ctools_safe_mul_size((size_t)nquantiles, 2 * sizeof(ST_double), &se_size) != 0) {
            return CBINSCATTER_ERR_MEMORY;
        }
        if (alloc_size > SIZE_MAX - se_size) {
            return CBINSCATTER_ERR_MEMORY;
        }
        alloc_size += se_size;
    }

    void *alloc_block = calloc(1, alloc_size);
    if (!alloc_block) {
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Partition the allocation */
    sum_y = (ST_double *)alloc_block;
    sum_x = sum_y + nquantiles;
    sum_w = sum_x + nquantiles;
    counts = (ST_int *)(sum_w + nquantiles);
    if (compute_se) {
        sum_y2 = (ST_double *)(counts + nquantiles);
        sum_x2 = sum_y2 + nquantiles;
    }

    /* Accumulate sums by bin - single pass through data */
    if (weights == NULL && !compute_se) {
        /* Fast path: unweighted, no SE */
        for (i = 0; i < N; i++) {
            b = bin_ids[i];
            if (b < 1 || b > nquantiles) continue;
            b--;
            sum_y[b] += y[i];
            sum_x[b] += x[i];
            counts[b]++;
        }
        /* sum_w = counts for unweighted */
        for (b = 0; b < nquantiles; b++) {
            sum_w[b] = (ST_double)counts[b];
        }
    } else if (weights == NULL) {
        /* Unweighted with SE */
        for (i = 0; i < N; i++) {
            b = bin_ids[i];
            if (b < 1 || b > nquantiles) continue;
            b--;
            yi = y[i];
            xi = x[i];
            sum_y[b] += yi;
            sum_x[b] += xi;
            sum_y2[b] += yi * yi;
            sum_x2[b] += xi * xi;
            counts[b]++;
        }
        for (b = 0; b < nquantiles; b++) {
            sum_w[b] = (ST_double)counts[b];
        }
    } else if (!compute_se) {
        /* Weighted, no SE */
        for (i = 0; i < N; i++) {
            b = bin_ids[i];
            if (b < 1 || b > nquantiles) continue;
            b--;
            w = weights[i];
            sum_y[b] += w * y[i];
            sum_x[b] += w * x[i];
            sum_w[b] += w;
            counts[b]++;
        }
    } else {
        /* Weighted with SE */
        for (i = 0; i < N; i++) {
            b = bin_ids[i];
            if (b < 1 || b > nquantiles) continue;
            b--;
            w = weights[i];
            yi = y[i];
            xi = x[i];
            sum_y[b] += w * yi;
            sum_x[b] += w * xi;
            sum_w[b] += w;
            sum_y2[b] += w * yi * yi;
            sum_x2[b] += w * xi * xi;
            counts[b]++;
        }
    }

    /* Count actual bins with data */
    for (b = 0; b < nquantiles; b++) {
        if (counts[b] > 0) actual_bins++;
    }

    /* Allocate result bins */
    result->num_bins = actual_bins;
    result->bins = (BinStats *)calloc(actual_bins, sizeof(BinStats));
    if (!result->bins) {
        free(alloc_block);
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Compute means and SE for each bin */
    ST_int bin_idx = 0;
    for (b = 0; b < nquantiles; b++) {
        if (counts[b] == 0) continue;

        result->bins[bin_idx].bin_id = b + 1;
        result->bins[bin_idx].n_obs = counts[b];
        result->bins[bin_idx].sum_weights = sum_w[b];

        /* Skip mean computation if sum of weights is zero */
        if (sum_w[b] <= 0.0) {
            result->bins[bin_idx].y_mean = 0.0;
            result->bins[bin_idx].x_mean = 0.0;
            result->bins[bin_idx].y_se = 0.0;
            result->bins[bin_idx].x_se = 0.0;
            bin_idx++;
            continue;
        }

        /* Compute weighted means */
        mean_y = sum_y[b] / sum_w[b];
        mean_x = sum_x[b] / sum_w[b];
        result->bins[bin_idx].y_mean = mean_y;
        result->bins[bin_idx].x_mean = mean_x;

        /* Compute SE if requested */
        if (compute_se && counts[b] > 1) {
            var_y = (sum_y2[b] / sum_w[b]) - (mean_y * mean_y);
            var_x = (sum_x2[b] / sum_w[b]) - (mean_x * mean_x);

            if (var_y > 0) {
                result->bins[bin_idx].y_se = sqrt(var_y * counts[b] / (counts[b] - 1) / counts[b]);
            }
            if (var_x > 0) {
                result->bins[bin_idx].x_se = sqrt(var_x * counts[b] / (counts[b] - 1) / counts[b]);
            }
        }

        bin_idx++;
    }

    free(alloc_block);
    return CBINSCATTER_OK;
}

/* ========================================================================
 * Single Group Bin Computation (Optimized)
 *
 * Uses quickselect O(N) for unweighted case, full sort O(N log N) for weighted.
 * ======================================================================== */

ST_retcode compute_bins_single_group(
    ST_double *y,
    ST_double *x,
    ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    ByGroupResult *result
) {
    ST_retcode rc = CBINSCATTER_OK;
    ST_int *sort_idx = NULL;
    ST_int *bin_ids = NULL;
    ST_int i, nq;

    nq = config->nquantiles;
    if (N < nq) nq = N;
    if (N < 2) return CBINSCATTER_ERR_FEW_OBS;

    /* Allocate sort indices and bin IDs */
    sort_idx = (ST_int *)malloc(N * sizeof(ST_int));
    bin_ids = (ST_int *)malloc(N * sizeof(ST_int));
    if (!sort_idx || !bin_ids) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Initialize sort indices (8x unrolled) */
    ST_int N8 = N - (N % 8);
    for (i = 0; i < N8; i += 8) {
        sort_idx[i]     = i;
        sort_idx[i + 1] = i + 1;
        sort_idx[i + 2] = i + 2;
        sort_idx[i + 3] = i + 3;
        sort_idx[i + 4] = i + 4;
        sort_idx[i + 5] = i + 5;
        sort_idx[i + 6] = i + 6;
        sort_idx[i + 7] = i + 7;
    }
    for (; i < N; i++) {
        sort_idx[i] = i;
    }

    if (weights == NULL && N >= HISTOGRAM_THRESHOLD) {
        /*
         * Unweighted, large N: Use histogram-based O(N) binning
         *
         * Strategy:
         * 1. Find min/max of x in one pass
         * 2. Build histogram of x values (HISTOGRAM_BUCKETS buckets)
         * 3. Use cumulative histogram to find bin boundaries
         * 4. Assign bins based on histogram bucket mapping
         *
         * This has excellent cache behavior - all sequential access.
         */
        ST_int *histogram = NULL;
        ST_int *cum_hist = NULL;
        ST_int *bucket_to_bin = NULL;

        histogram = (ST_int *)calloc(HISTOGRAM_BUCKETS, sizeof(ST_int));
        cum_hist = (ST_int *)malloc(HISTOGRAM_BUCKETS * sizeof(ST_int));
        bucket_to_bin = (ST_int *)malloc(HISTOGRAM_BUCKETS * sizeof(ST_int));

        if (!histogram || !cum_hist || !bucket_to_bin) {
            free(histogram);
            free(cum_hist);
            free(bucket_to_bin);
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        /* Pass 1: Find min/max - must skip missing values for initialization */
        ST_double x_min = 0, x_max = 0;
        int found_valid = 0;
        for (i = 0; i < N; i++) {
            if (!IS_MISSING(x[i])) {
                if (!found_valid) {
                    x_min = x[i];
                    x_max = x[i];
                    found_valid = 1;
                } else {
                    if (x[i] < x_min) x_min = x[i];
                    if (x[i] > x_max) x_max = x[i];
                }
            }
        }

        /* Handle all-missing case */
        if (!found_valid) {
            for (i = 0; i < N; i++) {
                bin_ids[i] = 0;
            }
            free(histogram);
            free(cum_hist);
            free(bucket_to_bin);
            goto compute_stats;
        }

        ST_double range = x_max - x_min;
        if (range <= 0.0) {
            /* All values identical - put everything in bin 1 */
            for (i = 0; i < N; i++) {
                bin_ids[i] = IS_MISSING(x[i]) ? 0 : 1;
            }
            free(histogram);
            free(cum_hist);
            free(bucket_to_bin);
            goto compute_stats;
        }

        ST_double scale = (HISTOGRAM_BUCKETS - 1) / range;

        /* Pass 2: Build histogram */
        for (i = 0; i < N; i++) {
            if (!IS_MISSING(x[i])) {
                ST_int bucket = (ST_int)((x[i] - x_min) * scale);
                if (bucket >= HISTOGRAM_BUCKETS) bucket = HISTOGRAM_BUCKETS - 1;
                if (bucket < 0) bucket = 0;
                histogram[bucket]++;
            }
        }

        /* Compute cumulative histogram and map buckets to bins */
        cum_hist[0] = histogram[0];
        for (i = 1; i < HISTOGRAM_BUCKETS; i++) {
            cum_hist[i] = cum_hist[i - 1] + histogram[i];
        }

        ST_int total_valid = cum_hist[HISTOGRAM_BUCKETS - 1];

        /* Map each histogram bucket to a quantile bin */
        for (i = 0; i < HISTOGRAM_BUCKETS; i++) {
            if (cum_hist[i] == 0) {
                bucket_to_bin[i] = 1;
            } else {
                ST_int bin = ((cum_hist[i] - 1) * nq) / total_valid + 1;
                if (bin > nq) bin = nq;
                if (bin < 1) bin = 1;
                bucket_to_bin[i] = bin;
            }
        }

        /* Pass 3: Assign bins using bucket lookup */
        for (i = 0; i < N; i++) {
            if (IS_MISSING(x[i])) {
                bin_ids[i] = 0;
            } else {
                ST_int bucket = (ST_int)((x[i] - x_min) * scale);
                if (bucket >= HISTOGRAM_BUCKETS) bucket = HISTOGRAM_BUCKETS - 1;
                if (bucket < 0) bucket = 0;
                bin_ids[i] = bucket_to_bin[bucket];
            }
        }

        free(histogram);
        free(cum_hist);
        free(bucket_to_bin);

    } else if (weights == NULL) {
        /*
         * Unweighted, small N: Full sort is fine
         * Assign bins directly from sorted order.
         */
        g_sort_values = x;
        qsort(sort_idx, N, sizeof(ST_int), compare_indices);

        for (i = 0; i < N; i++) {
            ST_int orig_idx = sort_idx[i];
            if (IS_MISSING(x[orig_idx])) {
                bin_ids[orig_idx] = 0;
            } else {
                ST_int bin = (ST_int)(((int64_t)i * nq) / N) + 1;
                if (bin > nq) bin = nq;
                bin_ids[orig_idx] = bin;
            }
        }

    } else if (N >= HISTOGRAM_THRESHOLD) {
        /*
         * Weighted, large N: Use histogram-based O(N) binning
         * Similar to unweighted but accumulate weights per bucket
         */
        ST_double *bucket_weights = NULL;
        ST_double *cum_weights = NULL;
        ST_int *bucket_to_bin = NULL;

        bucket_weights = (ST_double *)calloc(HISTOGRAM_BUCKETS, sizeof(ST_double));
        cum_weights = (ST_double *)malloc(HISTOGRAM_BUCKETS * sizeof(ST_double));
        bucket_to_bin = (ST_int *)malloc(HISTOGRAM_BUCKETS * sizeof(ST_int));

        if (!bucket_weights || !cum_weights || !bucket_to_bin) {
            free(bucket_weights);
            free(cum_weights);
            free(bucket_to_bin);
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        /* Pass 1: Find min/max - must skip missing values for initialization */
        ST_double x_min = 0, x_max = 0;
        int found_valid = 0;
        for (i = 0; i < N; i++) {
            if (!IS_MISSING(x[i])) {
                if (!found_valid) {
                    x_min = x[i];
                    x_max = x[i];
                    found_valid = 1;
                } else {
                    if (x[i] < x_min) x_min = x[i];
                    if (x[i] > x_max) x_max = x[i];
                }
            }
        }

        /* Handle all-missing case */
        if (!found_valid) {
            for (i = 0; i < N; i++) {
                bin_ids[i] = 0;
            }
            free(bucket_weights);
            free(cum_weights);
            free(bucket_to_bin);
            goto compute_stats;
        }

        ST_double range = x_max - x_min;
        if (range <= 0.0) {
            for (i = 0; i < N; i++) {
                bin_ids[i] = IS_MISSING(x[i]) ? 0 : 1;
            }
            free(bucket_weights);
            free(cum_weights);
            free(bucket_to_bin);
            goto compute_stats;
        }

        ST_double scale = (HISTOGRAM_BUCKETS - 1) / range;

        /* Pass 2: Accumulate weights per bucket */
        for (i = 0; i < N; i++) {
            if (!IS_MISSING(x[i])) {
                ST_int bucket = (ST_int)((x[i] - x_min) * scale);
                if (bucket >= HISTOGRAM_BUCKETS) bucket = HISTOGRAM_BUCKETS - 1;
                if (bucket < 0) bucket = 0;
                bucket_weights[bucket] += weights[i];
            }
        }

        /* Compute cumulative weights and map buckets to bins */
        cum_weights[0] = bucket_weights[0];
        for (i = 1; i < HISTOGRAM_BUCKETS; i++) {
            cum_weights[i] = cum_weights[i - 1] + bucket_weights[i];
        }

        ST_double total_weight = cum_weights[HISTOGRAM_BUCKETS - 1];

        /* Handle zero total weight - assign all to bin 1 */
        if (total_weight <= 0.0) {
            for (i = 0; i < N; i++) {
                bin_ids[i] = IS_MISSING(x[i]) ? 0 : 1;
            }
            free(bucket_weights);
            free(cum_weights);
            free(bucket_to_bin);
            goto compute_stats;
        }

        ST_double weight_per_bin = total_weight / nq;

        /* Map buckets to bins based on cumulative weight */
        for (i = 0; i < HISTOGRAM_BUCKETS; i++) {
            if (cum_weights[i] <= 0.0) {
                bucket_to_bin[i] = 1;
            } else {
                ST_int bin = (ST_int)((cum_weights[i] - 0.5 * bucket_weights[i]) / weight_per_bin) + 1;
                if (bin > nq) bin = nq;
                if (bin < 1) bin = 1;
                bucket_to_bin[i] = bin;
            }
        }

        /* Pass 3: Assign bins */
        for (i = 0; i < N; i++) {
            if (IS_MISSING(x[i])) {
                bin_ids[i] = 0;
            } else {
                ST_int bucket = (ST_int)((x[i] - x_min) * scale);
                if (bucket >= HISTOGRAM_BUCKETS) bucket = HISTOGRAM_BUCKETS - 1;
                if (bucket < 0) bucket = 0;
                bin_ids[i] = bucket_to_bin[bucket];
            }
        }

        free(bucket_weights);
        free(cum_weights);
        free(bucket_to_bin);

    } else {
        /*
         * Weighted, small N: Full sort for exact results
         */
        g_sort_values = x;
        qsort(sort_idx, N, sizeof(ST_int), compare_indices);

        ST_double total_weight = 0.0;
        ST_double *cum_weight = (ST_double *)malloc(N * sizeof(ST_double));
        if (!cum_weight) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        for (i = 0; i < N; i++) {
            ST_int orig_idx = sort_idx[i];
            if (!IS_MISSING(x[orig_idx])) {
                total_weight += weights[orig_idx];
            }
            cum_weight[i] = total_weight;
        }

        /* Handle zero total weight - assign all to bin 1 */
        if (total_weight <= 0.0) {
            for (i = 0; i < N; i++) {
                bin_ids[i] = IS_MISSING(x[i]) ? 0 : 1;
            }
            free(cum_weight);
            goto compute_stats;
        }

        ST_double weight_per_bin = total_weight / nq;
        for (i = 0; i < N; i++) {
            ST_int orig_idx = sort_idx[i];
            if (IS_MISSING(x[orig_idx])) {
                bin_ids[orig_idx] = 0;
            } else {
                ST_int bin = (ST_int)(cum_weight[i] / weight_per_bin) + 1;
                if (bin > nq) bin = nq;
                if (bin < 1) bin = 1;
                bin_ids[orig_idx] = bin;
            }
        }
        free(cum_weight);
    }

compute_stats:
    /* Compute bin statistics */
    rc = compute_bin_statistics(y, x, bin_ids, weights, N, nq,
                                config->compute_se, result);

cleanup:
    free(sort_idx);
    free(bin_ids);
    return rc;
}

/* ========================================================================
 * Discrete Mode Bin Computation
 * ======================================================================== */

ST_retcode compute_bins_discrete(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int compute_se,
    ByGroupResult *result
) {
    ST_int i, j, unique_count;
    ST_int *sort_idx = NULL;
    ST_double *unique_x = NULL;
    ST_int *bin_ids = NULL;
    ST_retcode rc = CBINSCATTER_OK;

    /* Allocate sort indices */
    sort_idx = (ST_int *)malloc(N * sizeof(ST_int));
    if (!sort_idx) return CBINSCATTER_ERR_MEMORY;

    for (i = 0; i < N; i++) {
        sort_idx[i] = i;
    }

    /* Argsort by x */
    g_sort_values = x;
    qsort(sort_idx, N, sizeof(ST_int), compare_indices);

    /* Count unique non-missing values */
    unique_count = 0;
    ST_double prev_val = STATA_MISSING;
    for (i = 0; i < N; i++) {
        ST_int idx = sort_idx[i];
        ST_double val = x[idx];
        if (IS_MISSING(val)) continue;
        if (unique_count == 0 || val != prev_val) {
            unique_count++;
            prev_val = val;
        }
    }

    if (unique_count < 1) {
        free(sort_idx);
        return CBINSCATTER_ERR_NOOBS;
    }

    /* Collect unique x values and assign bins in one pass */
    unique_x = (ST_double *)malloc(unique_count * sizeof(ST_double));
    bin_ids = (ST_int *)calloc(N, sizeof(ST_int));
    if (!unique_x || !bin_ids) {
        free(sort_idx);
        free(unique_x);
        free(bin_ids);
        return CBINSCATTER_ERR_MEMORY;
    }

    j = 0;
    prev_val = STATA_MISSING;
    for (i = 0; i < N; i++) {
        ST_int idx = sort_idx[i];
        ST_double val = x[idx];
        if (IS_MISSING(val)) {
            bin_ids[idx] = 0;
            continue;
        }
        if (j == 0 || val != prev_val) {
            unique_x[j] = val;
            j++;
            prev_val = val;
        }
        bin_ids[idx] = j;  /* 1-based bin ID */
    }

    /* Compute statistics */
    rc = compute_bin_statistics(y, x, bin_ids, weights, N, unique_count,
                                compute_se, result);

    free(sort_idx);
    free(unique_x);
    free(bin_ids);
    return rc;
}

/* ========================================================================
 * Legacy functions for API compatibility
 * ======================================================================== */

ST_retcode compute_quantile_cutpoints(
    const ST_double *x_sorted,
    ST_int N,
    ST_int nquantiles,
    const ST_double *w_sorted,
    ST_double *cutpoints
) {
    ST_int q, idx;
    ST_double p;

    if (w_sorted == NULL) {
        /* Unweighted: simple percentile computation */
        for (q = 0; q <= nquantiles; q++) {
            p = (ST_double)q / nquantiles;
            idx = (ST_int)(p * (N - 1));
            if (idx >= N) idx = N - 1;
            if (idx < 0) idx = 0;
            cutpoints[q] = x_sorted[idx];
        }
    } else {
        /* Weighted: compute total, then find cutpoints via cumulative weight */
        ST_double total_weight = 0.0;
        for (idx = 0; idx < N; idx++) {
            total_weight += w_sorted[idx];
        }

        cutpoints[0] = x_sorted[0];
        cutpoints[nquantiles] = x_sorted[N - 1];

        /* Use binary search on cumulative weights for O(N + nq*log(N)) */
        ST_double *cum_w = (ST_double *)malloc(N * sizeof(ST_double));
        if (cum_w) {
            cum_w[0] = w_sorted[0];
            for (idx = 1; idx < N; idx++) {
                cum_w[idx] = cum_w[idx - 1] + w_sorted[idx];
            }

            for (q = 1; q < nquantiles; q++) {
                ST_double target = ((ST_double)q / nquantiles) * total_weight;
                /* Binary search for target weight */
                ST_int lo = 0, hi = N - 1;
                while (lo < hi) {
                    ST_int mid = (lo + hi) / 2;
                    if (cum_w[mid] < target) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                cutpoints[q] = x_sorted[lo];
            }
            free(cum_w);
        } else {
            /* Fallback to O(nq * N) if allocation fails */
            for (q = 1; q < nquantiles; q++) {
                ST_double target_weight = ((ST_double)q / nquantiles) * total_weight;
                ST_double cum_weight = 0.0;
                for (idx = 0; idx < N; idx++) {
                    cum_weight += w_sorted[idx];
                    if (cum_weight >= target_weight) {
                        cutpoints[q] = x_sorted[idx];
                        break;
                    }
                }
                if (idx >= N) {
                    cutpoints[q] = x_sorted[N - 1];
                }
            }
        }
    }

    return CBINSCATTER_OK;
}

void assign_bins(
    const ST_double *x,
    ST_int N,
    const ST_double *cutpoints,
    ST_int nquantiles,
    ST_int *bin_ids
) {
    ST_int i;

    for (i = 0; i < N; i++) {
        if (IS_MISSING(x[i])) {
            bin_ids[i] = 0;
        } else {
            /* Binary search for bin */
            ST_int lo = 1, hi = nquantiles;
            ST_double xi = x[i];

            /* Quick bounds check */
            if (xi <= cutpoints[1]) {
                bin_ids[i] = 1;
                continue;
            }
            if (xi > cutpoints[nquantiles - 1]) {
                bin_ids[i] = nquantiles;
                continue;
            }

            /* Binary search */
            while (lo < hi) {
                ST_int mid = (lo + hi) / 2;
                if (xi > cutpoints[mid]) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            bin_ids[i] = lo;
        }
    }
}

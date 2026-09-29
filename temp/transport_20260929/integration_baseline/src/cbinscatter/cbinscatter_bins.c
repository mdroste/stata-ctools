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
#include "../ctools_select.h"
#include "../ctools_types.h"

/* Use SF_is_missing() from stplugin.h for missing value checks */

/* ========================================================================
 * Exact sort of x for quantile cutpoints
 * ======================================================================== */

static inline ST_double sortable_to_double(uint64_t key)
{
    uint64_t bits = (key >> 63) ? (key ^ UINT64_C(0x8000000000000000)) : ~key;
    ST_double d;
    memcpy(&d, &bits, sizeof(d));
    return d;
}

/*
 * Copy the non-missing x values (and their weights, when w != NULL) into
 * x_sorted/w_sorted in ascending order of x. LSD radix sort on sortable 64-bit
 * keys: exact and O(N) for any distribution (ties, outliers, skew), with no
 * shared comparator state. Returns the number of values copied, or -1 when
 * memory is short.
 */
static ST_int sort_nonmissing_x(const ST_double *x, const ST_double *w, ST_int N,
                                ST_double *x_sorted, ST_double *w_sorted)
{
    ST_int n = 0;
    for (ST_int i = 0; i < N; i++) {
        if (!SF_is_missing(x[i])) n++;
    }
    if (n == 0) return 0;

    uint64_t *keys = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *keys_tmp = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    ST_double *w_tmp = w ? (ST_double *)malloc((size_t)n * sizeof(ST_double)) : NULL;
    size_t (*counts)[256] = calloc(8, sizeof(*counts));
    if (!keys || !keys_tmp || (w && !w_tmp) || !counts) {
        free(keys); free(keys_tmp); free(w_tmp); free(counts);
        return -1;
    }

    ST_int j = 0;
    for (ST_int i = 0; i < N; i++) {
        if (SF_is_missing(x[i])) continue;
        uint64_t k = ctools_double_to_sortable(x[i]);
        keys[j] = k;
        if (w) w_sorted[j] = w[i];
        for (int p = 0; p < 8; p++) counts[p][(k >> (8 * p)) & 0xFF]++;
        j++;
    }

    uint64_t *src = keys, *dst = keys_tmp;
    ST_double *wsrc = w_sorted, *wdst = w_tmp;
    for (int p = 0; p < 8; p++) {
        size_t *c = counts[p];
        /* Skip a pass in which every key has the same byte */
        if (c[(src[0] >> (8 * p)) & 0xFF] == (size_t)n) continue;
        size_t offset = 0;
        for (int b = 0; b < 256; b++) {
            size_t t = c[b];
            c[b] = offset;
            offset += t;
        }
        for (ST_int i = 0; i < n; i++) {
            size_t pos = c[(src[i] >> (8 * p)) & 0xFF]++;
            dst[pos] = src[i];
            if (w) wdst[pos] = wsrc[i];
        }
        uint64_t *t = src; src = dst; dst = t;
        ST_double *tw = wsrc; wsrc = wdst; wdst = tw;
    }

    for (ST_int i = 0; i < n; i++) x_sorted[i] = sortable_to_double(src[i]);
    if (w && wsrc != w_sorted) memcpy(w_sorted, wsrc, (size_t)n * sizeof(ST_double));

    free(keys); free(keys_tmp); free(w_tmp); free(counts);
    return n;
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
    if (SF_is_missing(va) && SF_is_missing(vb)) return 0;
    if (SF_is_missing(va)) return 1;
    if (SF_is_missing(vb)) return -1;

    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

/* ========================================================================
 * Bin Statistics (optimized single-pass)
 * ======================================================================== */

static ST_retcode compute_bin_statistics(
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
 * Single Group Bin Computation
 * ======================================================================== */

ST_retcode compute_bins_single_group(
    ST_double *y,
    ST_double *x,
    ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    ByGroupResult *result,
    ST_int *bin_ids_out
) {
    ST_retcode rc = CBINSCATTER_OK;
    ST_int *bin_ids = bin_ids_out;
    ST_double *x_sorted = NULL, *w_sorted = NULL, *cutpoints = NULL;
    ST_int i, nq;

    nq = config->nquantiles;
    if (N < 2) return CBINSCATTER_ERR_FEW_OBS;

    if (!bin_ids) bin_ids = (ST_int *)malloc(N * sizeof(ST_int));
    x_sorted = (ST_double *)malloc(N * sizeof(ST_double));
    if (weights) w_sorted = (ST_double *)malloc(N * sizeof(ST_double));
    if (!bin_ids || !x_sorted || (weights && !w_sorted)) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /*
     * Quantile bins for every N, as binscatter's fastxtile and xtile form
     * them: the cutpoints are _pctile's (weighted) percentiles of x, and each
     * value goes to the first bin whose cutpoint is >= x. Tied values share a
     * bin, so heavily tied x gives fewer than nquantiles non-empty bins.
     */
    ST_int n_valid = sort_nonmissing_x(x, weights, N, x_sorted, w_sorted);
    if (n_valid < 0) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }
    if (n_valid == 0) {
        for (i = 0; i < N; i++) bin_ids[i] = 0;
        nq = 1;
    } else {
        if (n_valid < nq) nq = n_valid;
        cutpoints = (ST_double *)malloc((nq + 1) * sizeof(ST_double));
        if (!cutpoints) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }
        compute_quantile_cutpoints(x_sorted, n_valid, nq, w_sorted, cutpoints);
        assign_bins(x, N, cutpoints, nq, bin_ids);
    }

    /* Compute bin statistics */
    rc = compute_bin_statistics(y, x, bin_ids, weights, N, nq,
                                config->compute_se, result);

cleanup:
    free(x_sorted);
    free(w_sorted);
    free(cutpoints);
    if (bin_ids != bin_ids_out) free(bin_ids);
    return rc;
}

/*
 * Renumber bin ids to 1..result->num_bins, the positions of the non-empty
 * bins in result->bins. Tied cutpoints leave empty bins, so the raw ids can
 * have gaps; the regression adjustments index result->bins by id.
 */
ST_retcode densify_bin_ids(ST_int *bin_ids, ST_int N, const ByGroupResult *result)
{
    if (result->num_bins <= 0 || result->bins == NULL) {
        for (ST_int i = 0; i < N; i++) bin_ids[i] = 0;
        return CBINSCATTER_OK;
    }
    ST_int max_id = result->bins[result->num_bins - 1].bin_id;
    ST_int *map = (ST_int *)calloc((size_t)max_id + 1, sizeof(ST_int));
    if (!map) return CBINSCATTER_ERR_MEMORY;
    for (ST_int r = 0; r < result->num_bins; r++) {
        map[result->bins[r].bin_id] = r + 1;
    }
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i];
        bin_ids[i] = (b >= 1 && b <= max_id) ? map[b] : 0;
    }
    free(map);
    return CBINSCATTER_OK;
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
    ByGroupResult *result,
    ST_int *bin_ids_out
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
    ST_double prev_val = SV_missval;
    for (i = 0; i < N; i++) {
        ST_int idx = sort_idx[i];
        ST_double val = x[idx];
        if (SF_is_missing(val)) continue;
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
    prev_val = SV_missval;
    for (i = 0; i < N; i++) {
        ST_int idx = sort_idx[i];
        ST_double val = x[idx];
        if (SF_is_missing(val)) {
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

    /* Every unique value is a non-empty bin, so these ids are already dense */
    if (bin_ids_out) memcpy(bin_ids_out, bin_ids, (size_t)N * sizeof(ST_int));

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

    if (w_sorted == NULL) {
        /*
         * Unweighted: match Stata's _pctile/xtile algorithm
         *
         * Stata's xtile uses: cutpoint[q] = x[ceil(q * N / nquantiles)]
         * This ensures exactly (q * N / nquantiles) observations are <= cutpoint[q]
         */
        cutpoints[0] = x_sorted[0];  /* Minimum */
        cutpoints[nquantiles] = x_sorted[N - 1];  /* Maximum */

        for (q = 1; q < nquantiles; q++) {
            /* Stata's formula: ceil(q * N / nquantiles) gives 1-based index */
            idx = (ST_int)(((int64_t)q * N + nquantiles - 1) / nquantiles);  /* Ceiling division */
            if (idx >= N) idx = N - 1;
            if (idx < 1) idx = 1;
            cutpoints[q] = x_sorted[idx - 1];  /* Convert to 0-based */
        }
    } else {
        /* Weighted: compute total, then find cutpoints via cumulative weight */
        ST_double total_weight = 0.0;
        for (idx = 0; idx < N; idx++) {
            total_weight += w_sorted[idx];
        }

        cutpoints[0] = x_sorted[0];
        cutpoints[nquantiles] = x_sorted[N - 1];

        /* Use binary search on cumulative weights for O(N + nq*log(N)).
         * Kahan summation for numerical stability with unbalanced weights. */
        ST_double *cum_w = (ST_double *)malloc(N * sizeof(ST_double));
        if (cum_w) {
            ST_double kahan_sum = 0.0, kahan_c = 0.0;
            for (idx = 0; idx < N; idx++) {
                ST_double y_val = w_sorted[idx] - kahan_c;
                ST_double t = kahan_sum + y_val;
                kahan_c = (t - kahan_sum) - y_val;
                kahan_sum = t;
                cum_w[idx] = kahan_sum;
            }

            for (q = 1; q < nquantiles; q++) {
                /* W*q/nq is exact whenever it is an integer (e.g. fweights), so
                 * a cumulative weight equal to it is found as a tie, as in the
                 * unweighted rule; (q/nq)*W rounds and can move a cutpoint */
                ST_double target = total_weight * q / nquantiles;
                /* Binary search for target weight */
                ST_int lo = 0, hi = N - 1;
                while (lo < hi) {
                    ST_int mid = lo + (hi - lo) / 2;
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
                ST_double target_weight = total_weight * q / nquantiles;
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
        if (SF_is_missing(x[i])) {
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

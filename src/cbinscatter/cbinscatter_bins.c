/*
 * cbinscatter_bins.c
 *
 * Bin computation implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 *
 * Optimizations:
 * - Exact quantile cutpoints by selection around sampled splitters, without
 *   sorting all of x (small samples are still sorted)
 * - Bin assignment and bin sums in one block-parallel pass
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
 * Bin Statistics (block-parallel single pass)
 * ======================================================================== */

/* Number of the n sorted values s[] that are < v (branchless search) */
static inline ST_int count_lt(const ST_double *s, ST_int n, ST_double v)
{
    if (n <= 0) return 0;
    const ST_double *b = s;
    for (ST_int len = n; len > 1; ) {
        ST_int half = len / 2;
        b = (b[half] < v) ? b + half : b;
        len -= half;
    }
    return (ST_int)(b - s) + (*b < v);
}

/* Number of the n sorted values s[] that are <= v (branchless search) */
static inline ST_int count_le(const ST_double *s, ST_int n, ST_double v)
{
    if (n <= 0) return 0;
    const ST_double *b = s;
    for (ST_int len = n; len > 1; ) {
        ST_int half = len / 2;
        b = (b[half] <= v) ? b + half : b;
        len -= half;
    }
    return (ST_int)(b - s) + (*b <= v);
}

/* Bin of xi for cutpoints cut[0..nq]: the first bin whose upper cutpoint is
 * >= xi, as assign_bins() assigns it (0 = missing). The cutpoints are sorted,
 * so that bin is one more than the number of cut[1..nq-1] below xi. */
static inline ST_int bin_of(ST_double xi, const ST_double *cut, ST_int nq)
{
    if (SF_is_missing(xi)) return 0;
    return 1 + count_lt(cut + 1, nq - 1, xi);
}

/*
 * Bin means (and, with compute_se, standard errors) of y and x. A row's bin
 * is bin_ids[i] when cutpoints is NULL; otherwise it is assigned from the
 * cutpoints (nquantiles + 1 values) and stored to bin_ids when that is not
 * NULL. Sums are formed per block of rows and added in block order.
 */
static ST_retcode compute_bin_statistics(
    const ST_double *y,
    const ST_double *x,
    ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int nquantiles,
    const ST_double *cutpoints,
    ST_int compute_se,
    ByGroupResult *result
) {
    ST_int b;
    ST_double mean_y, mean_x, var_y, var_x;
    ST_int actual_bins = 0;

    /* Validate nquantiles to prevent overflow */
    if (nquantiles <= 0 || (!cutpoints && !bin_ids)) {
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Per block: sum_y, sum_x, sum_w [, sum_y2, sum_x2] and counts, then the
     * same for the totals in one more slot (index num_blocks) */
    const ST_int num_blocks = cbs_num_blocks(N);
    const size_t nfields = compute_se ? 5 : 3;
    size_t per_block, n_dbl, n_int;
    if (ctools_safe_mul_size(nfields, (size_t)nquantiles, &per_block) != 0 ||
        ctools_safe_mul_size((size_t)num_blocks + 1, per_block, &n_dbl) != 0 ||
        ctools_safe_mul_size((size_t)num_blocks + 1, (size_t)nquantiles, &n_int) != 0) {
        return CBINSCATTER_ERR_MEMORY;
    }
    ST_double *sums = (ST_double *)calloc(n_dbl, sizeof(ST_double));
    ST_int *block_counts = (ST_int *)calloc(n_int, sizeof(ST_int));
    if (!sums || !block_counts) {
        free(sums); free(block_counts);
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Accumulate sums by bin - one pass through the data, blocks in parallel */
    ST_int blk;
    #pragma omp parallel for schedule(dynamic) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        ST_double *sum_y = sums + (size_t)blk * per_block;
        ST_double *sum_x = sum_y + nquantiles;
        ST_double *sum_w = sum_x + nquantiles;
        ST_double *sum_y2 = compute_se ? sum_w + nquantiles : NULL;
        ST_double *sum_x2 = compute_se ? sum_y2 + nquantiles : NULL;
        ST_int *counts = block_counts + (size_t)blk * nquantiles;
        for (ST_int i = i0; i < i1; i++) {
            ST_int bi;
            if (cutpoints) {
                bi = bin_of(x[i], cutpoints, nquantiles);
                if (bin_ids) bin_ids[i] = bi;
            } else {
                bi = bin_ids[i];
            }
            if (bi < 1 || bi > nquantiles) continue;
            bi--;
            /* Unweighted rows have weight 1: w * y is y and sum_w counts them */
            ST_double w = weights ? weights[i] : 1.0;
            ST_double yi = y[i], xi = x[i];
            sum_y[bi] += w * yi;
            sum_x[bi] += w * xi;
            sum_w[bi] += w;
            if (compute_se) {
                sum_y2[bi] += w * yi * yi;
                sum_x2[bi] += w * xi * xi;
            }
            counts[bi]++;
        }
    }

    /* Totals, adding the blocks in order */
    ST_double *sum_y = sums + (size_t)num_blocks * per_block;
    ST_double *sum_x = sum_y + nquantiles;
    ST_double *sum_w = sum_x + nquantiles;
    ST_double *sum_y2 = compute_se ? sum_w + nquantiles : NULL;
    ST_double *sum_x2 = compute_se ? sum_y2 + nquantiles : NULL;
    ST_int *counts = block_counts + (size_t)num_blocks * nquantiles;
    for (blk = 0; blk < num_blocks; blk++) {
        const ST_double *src = sums + (size_t)blk * per_block;
        const ST_int *csrc = block_counts + (size_t)blk * nquantiles;
        for (size_t f = 0; f < per_block; f++) sum_y[f] += src[f];
        for (b = 0; b < nquantiles; b++) counts[b] += csrc[b];
    }

    /* Count actual bins with data */
    for (b = 0; b < nquantiles; b++) {
        if (counts[b] > 0) actual_bins++;
    }

    /* Allocate result bins */
    result->num_bins = actual_bins;
    result->bins = (BinStats *)calloc(actual_bins, sizeof(BinStats));
    if (!result->bins) {
        free(sums); free(block_counts);
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

    free(sums);
    free(block_counts);
    return CBINSCATTER_OK;
}

/* ========================================================================
 * Quantile Cutpoints
 * ======================================================================== */

/* Below this many non-missing values, x is simply sorted */
#define CBS_SORT_MAX_N 131072
/* Values sampled to place the splitters */
#define CBS_SAMPLE_N 65536

/* Interval of v among the sorted distinct splitters B[0..nb-1]: 2k for the
 * open interval (B[k-1], B[k]) (unbounded at the ends) and 2k - 1 for the
 * value B[k-1] itself. Intervals 0..2nb are in increasing order of value. */
static inline ST_int interval_of(const ST_double *B, ST_int nb, ST_double v)
{
    ST_int k = count_le(B, nb, v);
    return (k > 0 && v == B[k - 1]) ? 2 * k - 1 : 2 * k;
}

static int compare_doubles(const void *a, const void *b)
{
    ST_double va = *(const ST_double *)a, vb = *(const ST_double *)b;
    return (va > vb) - (va < vb);
}

/*
 * The cutpoints compute_quantile_cutpoints() takes from the sorted
 * non-missing x (weighted when w is not NULL), without sorting all of x.
 * A sorted sample of x places two splitters around each cutpoint's rank (or
 * cumulative weight). One block-parallel pass counts the values in each
 * interval between and at the splitters, which fixes the interval that holds
 * each cutpoint; a second copies out only those intervals, and a selection
 * (weighted: a sort) within each gives the exact cutpoint. Splitters that
 * miss only enlarge the intervals copied. n is the number of non-missing x
 * (nq <= n), xmin and xmax their extremes. Returns 0, or -1 when memory is
 * short.
 */
static int cutpoints_by_selection(const ST_double *x, const ST_double *w, ST_int N, ST_int n,
                                  ST_double xmin, ST_double xmax, ST_int nq, ST_double *cut)
{
    int status = -1;
    const ST_int num_blocks = cbs_num_blocks(N);
    ST_double *sample = NULL, *sample_w = NULL, *sorted = NULL, *sorted_w = NULL;
    ST_double *B = NULL, *block_w = NULL, *tot_w = NULL, *cum_w = NULL;
    ST_double *buf = NULL, *buf_w = NULL, *targets_w = NULL;
    ST_int *block_n = NULL, *region_of = NULL;
    int64_t *tot = NULL, *start = NULL, *ranks = NULL, *base = NULL;
    ST_int *needed = NULL;

    cut[0] = xmin;
    cut[nq] = xmax;
    if (nq < 2) return 0;

    /* Sample at evenly spaced rows */
    ST_int s_max = CBS_SAMPLE_N, s = 0;
    sample = (ST_double *)malloc((size_t)s_max * sizeof(ST_double));
    sorted = (ST_double *)malloc((size_t)s_max * sizeof(ST_double));
    if (w) {
        sample_w = (ST_double *)malloc((size_t)s_max * sizeof(ST_double));
        sorted_w = (ST_double *)malloc((size_t)s_max * sizeof(ST_double));
    }
    if (!sample || !sorted || (w && (!sample_w || !sorted_w))) goto done;
    for (ST_int k = 0; k < s_max; k++) {
        ST_int i = (ST_int)(((2 * (int64_t)k + 1) * N) / (2 * (int64_t)s_max));
        if (SF_is_missing(x[i])) continue;
        sample[s] = x[i];
        if (w) sample_w[s] = w[i];
        s++;
    }
    if (s > 0 && sort_nonmissing_x(sample, sample_w, s, sorted, sorted_w) != s) goto done;

    /* Two splitters per cutpoint, about 3.5 standard errors of the sample
     * quantile either side of its estimate */
    B = (ST_double *)malloc((size_t)2 * nq * sizeof(ST_double));
    ranks = (int64_t *)malloc((size_t)nq * sizeof(int64_t));
    targets_w = (ST_double *)malloc((size_t)nq * sizeof(ST_double));
    if (!B || !ranks || !targets_w) goto done;
    ST_double sample_total = 0.0, sample_cum = 0.0;
    if (w) for (ST_int k = 0; k < s; k++) sample_total += sorted_w[k];
    ST_int nb = 0, pos = 0;
    for (ST_int q = 1; q < nq && s > 0; q++) {
        ST_double p = (ST_double)q / nq, e;
        if (w) {
            ST_double target = sample_total * q / nq;
            while (pos < s - 1 && sample_cum + sorted_w[pos] < target) sample_cum += sorted_w[pos++];
            e = pos;
        } else {
            e = p * s - 0.5;
        }
        ST_double delta = (w ? 5.0 : 3.5) * sqrt(s * p * (1.0 - p)) + 2.0;
        int64_t lo = (int64_t)floor(e - delta), hi = (int64_t)ceil(e + delta);
        if (lo >= 0) B[nb++] = sorted[lo];
        if (hi < s) B[nb++] = sorted[hi];
    }
    qsort(B, (size_t)nb, sizeof(ST_double), compare_doubles);
    {
        ST_int u = 0;
        for (ST_int k = 0; k < nb; k++) if (u == 0 || B[k] != B[u - 1]) B[u++] = B[k];
        nb = u;
    }
    const ST_int nr = 2 * nb + 1;

    /* Pass 1: values (and weights) per interval, per block */
    size_t cells;
    if (ctools_safe_mul_size((size_t)num_blocks, (size_t)nr, &cells) != 0) goto done;
    block_n = (ST_int *)calloc(cells, sizeof(ST_int));
    if (w) block_w = (ST_double *)calloc(cells, sizeof(ST_double));
    if (!block_n || (w && !block_w)) goto done;
    ST_int blk;
    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        ST_int *cn = block_n + (size_t)blk * nr;
        ST_double *cw = w ? block_w + (size_t)blk * nr : NULL;
        for (ST_int i = i0; i < i1; i++) {
            if (SF_is_missing(x[i])) continue;
            ST_int r = interval_of(B, nb, x[i]);
            cn[r]++;
            if (cw) cw[r] += w[i];
        }
    }

    /* Interval totals (blocks in order) and their starting ranks/weights */
    tot = (int64_t *)calloc((size_t)nr, sizeof(int64_t));
    start = (int64_t *)calloc((size_t)nr + 1, sizeof(int64_t));
    needed = (ST_int *)malloc((size_t)nr * sizeof(ST_int));
    region_of = (ST_int *)malloc((size_t)nq * sizeof(ST_int));
    if (w) {
        tot_w = (ST_double *)calloc((size_t)nr, sizeof(ST_double));
        cum_w = (ST_double *)calloc((size_t)nr + 1, sizeof(ST_double));
    }
    if (!tot || !start || !needed || !region_of || (w && (!tot_w || !cum_w))) goto done;
    for (blk = 0; blk < num_blocks; blk++) {
        for (ST_int r = 0; r < nr; r++) {
            tot[r] += block_n[(size_t)blk * nr + r];
            if (w) tot_w[r] += block_w[(size_t)blk * nr + r];
        }
    }
    for (ST_int r = 0; r < nr; r++) {
        start[r + 1] = start[r] + tot[r];
        if (w) cum_w[r + 1] = cum_w[r] + tot_w[r];
        needed[r] = -1;
    }

    /* The interval of each cutpoint. A value interval (odd r) is the
     * cutpoint itself; an open interval must be searched. */
    ST_double total_weight = w ? cum_w[nr] : 0.0;
    ST_int r = 0;
    for (ST_int q = 1; q < nq; q++) {
        region_of[q] = -1;
        if (w) {
            /* First value whose cumulative weight reaches W * q / nq */
            ST_double target = total_weight * q / nq;
            targets_w[q] = target;
            while (r < nr && !(tot[r] > 0 && cum_w[r + 1] >= target)) r++;
            if (r == nr) {
                cut[q] = xmax;
                r = nr - 1;
                continue;
            }
        } else {
            /* Rank ceil(q * n / nq), 1-based, as in compute_quantile_cutpoints */
            int64_t idx = ((int64_t)q * n + nq - 1) / nq;
            if (idx >= n) idx = n - 1;
            if (idx < 1) idx = 1;
            ranks[q] = idx - 1;
            while (r < nr && start[r + 1] <= ranks[q]) r++;
        }
        if (r % 2 == 1) {
            cut[q] = B[(r - 1) / 2];
        } else {
            region_of[q] = r;
            needed[r] = 0;
        }
    }

    /* Pass 2: copy out the needed intervals. base[r] is an interval's offset
     * in buf, block_n each block's write position within the interval. */
    base = (int64_t *)calloc((size_t)nr + 1, sizeof(int64_t));
    if (!base) goto done;
    int64_t buf_n = 0;
    for (r = 0; r < nr; r++) {
        base[r] = buf_n;
        if (needed[r] < 0) continue;
        int64_t run = buf_n;
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int c = block_n[(size_t)blk * nr + r];
            block_n[(size_t)blk * nr + r] = (ST_int)(run - buf_n);
            run += c;
        }
        buf_n = run;
    }
    if (buf_n > 0) {
        buf = (ST_double *)malloc((size_t)buf_n * sizeof(ST_double));
        if (w) buf_w = (ST_double *)malloc((size_t)buf_n * sizeof(ST_double));
        if (!buf || (w && !buf_w)) goto done;
        #pragma omp parallel for schedule(static) if (num_blocks > 1)
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int i0, i1;
            cbs_block_range(blk, N, &i0, &i1);
            ST_int *off = block_n + (size_t)blk * nr;
            for (ST_int i = i0; i < i1; i++) {
                if (SF_is_missing(x[i])) continue;
                ST_int ri = interval_of(B, nb, x[i]);
                if (needed[ri] < 0) continue;
                int64_t at = base[ri] + off[ri]++;
                buf[at] = x[i];
                if (w) buf_w[at] = w[i];
            }
        }
    }

    /* Select each cutpoint within its interval; intervals in parallel */
    int failed = 0;
    #pragma omp parallel for schedule(dynamic) reduction(| : failed) if (buf_n > CBS_SORT_MAX_N)
    for (ST_int rr = 0; rr < nr; rr++) {
        if (needed[rr] < 0) continue;
        ST_double *v = buf + base[rr];
        const int64_t len = tot[rr];
        if (!w) {
            /* Ascending ranks: each selection leaves larger values after it */
            int64_t lo = 0;
            for (ST_int q = 1; q < nq; q++) {
                if (region_of[q] != rr) continue;
                int64_t k = ranks[q] - start[rr];
                cut[q] = ctools_quickselect_double(v + lo, (size_t)(len - lo), (size_t)(k - lo));
                lo = k + 1;
            }
        } else {
            /* Sort the interval by x and accumulate its weights after those
             * of all smaller values */
            ST_double *xs = (ST_double *)malloc((size_t)len * sizeof(ST_double));
            ST_double *ws = (ST_double *)malloc((size_t)len * sizeof(ST_double));
            if (!xs || !ws || sort_nonmissing_x(v, buf_w + base[rr], (ST_int)len, xs, ws) != len) {
                failed |= 1;
            } else {
                /* Kahan-accumulated cumulative weight through xs[added - 1];
                 * ascending targets resume where the previous one stopped */
                ST_double kahan_sum = cum_w[rr], kahan_c = 0.0;
                int64_t added = 0;
                for (ST_int q = 1; q < nq; q++) {
                    if (region_of[q] != rr) continue;
                    while (added < len && (added == 0 || kahan_sum < targets_w[q])) {
                        ST_double y_val = ws[added] - kahan_c;
                        ST_double t = kahan_sum + y_val;
                        kahan_c = (t - kahan_sum) - y_val;
                        kahan_sum = t;
                        added++;
                    }
                    cut[q] = xs[added - 1];
                }
            }
            free(xs);
            free(ws);
        }
    }
    status = failed ? -1 : 0;

done:
    free(sample); free(sample_w); free(sorted); free(sorted_w);
    free(B); free(block_n); free(block_w); free(tot); free(start); free(tot_w); free(cum_w);
    free(needed); free(region_of); free(ranks); free(targets_w); free(base);
    free(buf); free(buf_w);
    return status;
}

/*
 * Quantile cutpoints of the non-missing x for nq bins: cut[0..*nq_out],
 * where *nq_out = min(nq, number of non-missing x), and 0 when every x is
 * missing. Returns CBINSCATTER_OK or CBINSCATTER_ERR_MEMORY.
 */
static ST_retcode quantile_cutpoints(const ST_double *x, const ST_double *weights, ST_int N,
                                     ST_int nq, ST_double *cut, ST_int *nq_out)
{
    /* Non-missing count and extremes, blocks in parallel */
    const ST_int num_blocks = cbs_num_blocks(N);
    ST_int n = 0;
    ST_double xmin = HUGE_VAL, xmax = -HUGE_VAL;
    ST_int blk;
    #pragma omp parallel for schedule(static) reduction(+ : n) reduction(min : xmin) \
        reduction(max : xmax) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        for (ST_int i = i0; i < i1; i++) {
            ST_double v = x[i];
            if (SF_is_missing(v)) continue;
            n++;
            if (v < xmin) xmin = v;
            if (v > xmax) xmax = v;
        }
    }
    *nq_out = 0;
    if (n == 0) return CBINSCATTER_OK;
    if (n < nq) nq = n;
    *nq_out = nq;

    if (n > CBS_SORT_MAX_N) {
        return cutpoints_by_selection(x, weights, N, n, xmin, xmax, nq, cut) == 0
                   ? CBINSCATTER_OK : CBINSCATTER_ERR_MEMORY;
    }

    ST_double *x_sorted = (ST_double *)malloc((size_t)n * sizeof(ST_double));
    ST_double *w_sorted = weights ? (ST_double *)malloc((size_t)n * sizeof(ST_double)) : NULL;
    ST_retcode rc = CBINSCATTER_OK;
    if (!x_sorted || (weights && !w_sorted) ||
        sort_nonmissing_x(x, weights, N, x_sorted, w_sorted) != n) {
        rc = CBINSCATTER_ERR_MEMORY;
    } else {
        compute_quantile_cutpoints(x_sorted, n, nq, w_sorted, cut);
    }
    free(x_sorted);
    free(w_sorted);
    return rc;
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
    const ST_double *cutpoints_in,
    ST_int nq_in,
    ByGroupResult *result,
    ST_int *bin_ids_out
) {
    ST_retcode rc;
    ST_double *cutpoints = NULL;
    ST_int nq = config->nquantiles;

    if (N < 1) return CBINSCATTER_ERR_FEW_OBS;

    if (cutpoints_in) {
        return compute_bin_statistics(y, x, bin_ids_out, weights, N, nq_in, cutpoints_in,
                                      config->compute_se, result);
    }

    /*
     * Quantile bins for every N, as binscatter's fastxtile and xtile form
     * them: the cutpoints are _pctile's (weighted) percentiles of x, and each
     * value goes to the first bin whose cutpoint is >= x. Tied values share a
     * bin, so heavily tied x gives fewer than nquantiles non-empty bins.
     */
    cutpoints = (ST_double *)malloc(((size_t)nq + 1) * sizeof(ST_double));
    if (!cutpoints) return CBINSCATTER_ERR_MEMORY;
    rc = quantile_cutpoints(x, weights, N, nq, cutpoints, &nq);
    if (rc == CBINSCATTER_OK && nq == 0) {
        /* Every x is missing: no bins */
        ST_int *zeros = (ST_int *)calloc((size_t)N, sizeof(ST_int));
        rc = zeros ? compute_bin_statistics(y, x, zeros, weights, N, 1, NULL,
                                            config->compute_se, result)
                   : CBINSCATTER_ERR_MEMORY;
        if (zeros && bin_ids_out) memcpy(bin_ids_out, zeros, (size_t)N * sizeof(ST_int));
        free(zeros);
    } else if (rc == CBINSCATTER_OK) {
        rc = compute_bin_statistics(y, x, bin_ids_out, weights, N, nq, cutpoints,
                                    config->compute_se, result);
    }
    free(cutpoints);
    return rc;
}

ST_retcode compute_pooled_cutpoints(
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int nquantiles,
    ST_double **cutpoints_out,
    ST_int *nq_out
) {
    *cutpoints_out = NULL;
    *nq_out = nquantiles;
    ST_double *cutpoints = (ST_double *)malloc(((size_t)nquantiles + 1) * sizeof(ST_double));
    if (!cutpoints) return CBINSCATTER_ERR_MEMORY;
    ST_int nq = 0;
    ST_retcode rc = quantile_cutpoints(x, weights, N, nquantiles, cutpoints, &nq);
    if (rc == CBINSCATTER_OK && nq > 0) {
        *cutpoints_out = cutpoints;
        *nq_out = nq;
    } else {
        free(cutpoints);
    }
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
    rc = compute_bin_statistics(y, x, bin_ids, weights, N, unique_count, NULL,
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

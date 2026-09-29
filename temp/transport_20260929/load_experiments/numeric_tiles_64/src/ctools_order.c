#include <stdlib.h>
#include <string.h>
#include "ctools_order.h"
#include "ctools_config.h"

int ctools_compare_rows(const stata_data *a, size_t i, const int *ak,
                       const stata_data *b, size_t j, const int *bk, size_t nk)
{
    for (size_t k = 0; k < nk; k++) {
        const stata_variable *x = a->vars + ak[k], *y = b->vars + bk[k];
        int c;
        if (x->type == STATA_TYPE_STRING) c = strcmp(x->data.str[i], y->data.str[j]);
        else {
            double u = x->data.dbl[i], v = y->data.dbl[j];
            c = (u > v) - (u < v);
        }
        if (c) return c;
    }
    return 0;
}

size_t ctools_group_starts(const stata_data *d, const perm_idx_t *order,
                           const int *keys, size_t nkeys, size_t *starts)
{
    size_t n = d->nobs, ngroups = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == 0 || ctools_compare_rows(d, order ? order[i - 1] : i - 1, keys,
                                          d, order ? order[i] : i, keys, nkeys))
            starts[ngroups++] = i;
    }
    starts[ngroups] = n;
    return ngroups;
}


static stata_retcode numeric_order(stata_data *d, const int *keys, size_t nk,
                                  perm_idx_t *tmp)
{
    size_t n = d->nobs;
    int nt = ctools_get_max_threads();
    if ((size_t)nt > (n + 16383) / 16384) nt = (int)((n + 16383) / 16384);
    uint64_t *bits = malloc(n * sizeof(*bits));
    size_t *hist = malloc((size_t)nt * 256 * sizeof(*hist));
    if (!bits || !hist) { free(bits); free(hist); return STATA_ERR_MEMORY; }
    perm_idx_t *src = d->sort_order, *dst = tmp;
    /* LSD over both columns and bytes gives a stable lexicographic order. */
    for (size_t k = nk; k-- > 0;) {
        const double *x = d->vars[keys[k]].data.dbl;
        uint64_t first = ctools_double_to_sortable(x[0]), varying = 0;
        #pragma omp parallel for num_threads(nt) reduction(|:varying)
        for (size_t i = 0; i < n; i++) {
            bits[i] = ctools_double_to_sortable(x[i]);
            varying |= bits[i] ^ first;
        }
        for (unsigned shift = 0; shift < 64; shift += 8) {
            if (!((varying >> shift) & 255)) continue;
            #pragma omp parallel num_threads(nt)
            {
                int tid = 0, team = 1;
                #ifdef _OPENMP
                tid = omp_get_thread_num(); team = omp_get_num_threads();
                #endif
                size_t begin = n * (size_t)tid / team, end = n * (size_t)(tid+1) / team;
                size_t *h = hist + (size_t)tid * 256;
                memset(h, 0, 256 * sizeof(*h));
                for (size_t i = begin; i < end; i++) h[(bits[src[i]] >> shift) & 255]++;
                #pragma omp barrier
                #pragma omp single
                {
                    size_t offset = 0;
                    for (size_t b = 0; b < 256; b++) {
                        for (int t = 0; t < team; t++) {
                            size_t *v = hist + (size_t)t * 256 + b, count = *v;
                            *v = offset; offset += count;
                        }
                    }
                }
                for (size_t i = begin; i < end; i++) {
                    perm_idx_t row = src[i];
                    dst[h[(bits[row] >> shift) & 255]++] = row;
                }
            }
            perm_idx_t *swap = src; src = dst; dst = swap;
        }
    }
    if (src != d->sort_order) memcpy(d->sort_order, src, n * sizeof(*src));
    free(bits); free(hist);
    return STATA_OK;
}

/* Locate a diagonal in a stable merge. Equal keys from the left run precede
 * keys from the right run. This lets the final, largest merge use every worker. */
static size_t merge_left(const stata_data *d, const int *keys, size_t nk,
                         const perm_idx_t *a, size_t na, const perm_idx_t *b,
                         size_t nb, size_t diagonal)
{
    size_t lo = diagonal > nb ? diagonal - nb : 0;
    size_t hi = diagonal < na ? diagonal : na;
    while (lo < hi) {
        size_t i = lo + (hi-lo)/2, j = diagonal-i;
        if (j && i < na && ctools_compare_rows(d, a[i], keys, d, b[j-1], keys, nk) <= 0)
            lo = i+1;
        else hi = i;
    }
    return lo;
}

stata_retcode ctools_order_stable(stata_data *d, const int *keys, size_t nk)
{
    size_t n = d->nobs;
    if (n < 2 || !nk) return STATA_OK;
    /* Sorted inputs are common for panel data and require no scratch space. */
    size_t i = 1;
    for (; i < n; i++)
        if (ctools_compare_rows(d, d->sort_order[i-1], keys, d, d->sort_order[i], keys, nk) > 0) break;
    if (i == n) return STATA_OK;
    perm_idx_t *tmp = malloc(n * sizeof(*tmp));
    if (!tmp) return STATA_ERR_MEMORY;
    int numeric = n >= 32768;
    for (size_t k = 0; k < nk; k++)
        if (d->vars[keys[k]].type != STATA_TYPE_DOUBLE) numeric = 0;
    if (numeric) {
        stata_retcode rc = numeric_order(d, keys, nk, tmp);
        free(tmp);
        return rc;
    }
    perm_idx_t *src = d->sort_order, *dst = tmp;
    for (size_t width = 1; width < n; width *= 2) {
        size_t stride = width * 2, tile = stride < 8192 ? stride : 8192;
        #pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if(n > 32768)
        for (size_t start = 0; start < n; start += tile) {
            size_t base = start / stride * stride;
            size_t mid = base + width < n ? base + width : n;
            size_t end = base + stride < n ? base + stride : n;
            size_t stop = start + tile < end ? start + tile : end;
            size_t l = base + merge_left(d, keys, nk, src+base, mid-base, src+mid, end-mid, start-base);
            size_t r = mid + start - base - (l-base), out = start;
            while (out < stop && l < mid && r < end) {
                if (ctools_compare_rows(d, src[l], keys, d, src[r], keys, nk) <= 0)
                    dst[out++] = src[l++];
                else dst[out++] = src[r++];
            }
            while (out < stop && l < mid) dst[out++] = src[l++];
            while (out < stop && r < end) dst[out++] = src[r++];
        }
        perm_idx_t *swap = src; src = dst; dst = swap;
    }
    if (src != d->sort_order) memcpy(d->sort_order, src, n * sizeof(*src));
    free(tmp);
    return STATA_OK;
}

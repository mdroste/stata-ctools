/*
 * creghdfe_core.c
 *
 * The regression's columns from the loaded variables, then the collinearity
 * rules and OLS on the partialled data (see creghdfe_core.h). Moved from
 * do_full_regression's STEPs 5 to 8 without changing the arithmetic.
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_core.h"

#include "../ctools_ols.h"
#include "../ctools_types.h"
#include "../ctools_unroll.h"
#include <stdint.h>

#if CTOOLS_UNROLL_K != 8
#error "creghdfe_columns_build reproduces the eight chains of ctools_dot_unrolled"
#endif

/* ========================================================================
 * The regression's columns from the loaded variables (creghdfe_core.h)
 * ======================================================================== */

#if defined(__x86_64__) || defined(_M_X64)
#include <emmintrin.h>
#define CREGHDFE_STREAM_STORES 1
#endif

/* Rows read from a loaded column at a time. A multiple of eight, so that a
 * block boundary never splits the chains of the sum of squares. */
#define CREGHDFE_COLUMN_BLOCK 1024

/* The matrix is written with stores that bypass the cache once it is larger
 * than this: an ordinary store reads every line it is about to overwrite.
 * Smaller matrices stay in the cache for the projection that follows. */
#define CREGHDFE_STREAM_BYTES ((size_t)32 << 20)

/* Column k on the kept rows i0 .. i0 + n - 1 (n <= CREGHDFE_COLUMN_BLOCK):
 * the loaded column itself where those are consecutive loaded rows, buf
 * (filled here) otherwise. */
static inline const ST_double *column_values(const creghdfe_columns *C, ST_int k, ST_int i0,
                                             ST_int n, ST_double *RESTRICT buf)
{
    const ST_double *RESTRICT src = C->cols[k];
    const ST_double level = C->level ? C->level[k] : -1.0;
    if (level >= 0)
    {
        if (C->rows)
        {
            const ST_int *RESTRICT rows = C->rows + i0;
            for (ST_int i = 0; i < n; i++)
                buf[i] = src[rows[i]] == level ? 1.0 : 0.0;
        }
        else
        {
            src += i0;
            for (ST_int i = 0; i < n; i++)
                buf[i] = src[i] == level ? 1.0 : 0.0;
        }
        return buf;
    }
    if (!C->rows)
        return src + i0;
    const ST_int *RESTRICT rows = C->rows + i0;
    for (ST_int i = 0; i < n; i++)
        buf[i] = src[rows[i]];
    return buf;
}

/* Sum of weight * value in row order */
static ST_double column_sum(const creghdfe_columns *C, ST_int k)
{
    const ST_int N = C->N;
    const ST_double *RESTRICT w = C->weights;
    ST_double buf[CREGHDFE_COLUMN_BLOCK];
    ST_double total = 0.0;
    for (ST_int i0 = 0; i0 < N; i0 += CREGHDFE_COLUMN_BLOCK)
    {
        const ST_int n = N - i0 < CREGHDFE_COLUMN_BLOCK ? N - i0 : CREGHDFE_COLUMN_BLOCK;
        const ST_double *RESTRICT x = column_values(C, k, i0, n, buf);
        if (w)
        {
            const ST_double *RESTRICT wi = w + i0;
            for (ST_int i = 0; i < n; i++)
                total += wi[i] * x[i];
        }
        else
        {
            for (ST_int i = 0; i < n; i++)
                total += x[i];
        }
    }
    return total;
}

void creghdfe_column_sums(const creghdfe_columns *C, ST_double *sums, ST_double *sum_w)
{
    const ST_int N = C->N, K = C->K;
    /* with weights one more task adds them up */
    const ST_int tasks = K + (C->weights && sum_w ? 1 : 0);
    ST_int k;
#pragma omp parallel for schedule(static) if (tasks > 1)
    for (k = 0; k < tasks; k++)
    {
        if (k == K)
        {
            const ST_double *RESTRICT w = C->weights;
            ST_double total = 0.0;
            for (ST_int i = 0; i < N; i++)
                total += w[i];
            *sum_w = total;
            continue;
        }
        sums[k] = column_sum(C, k);
    }
}

/* Compensated (weighted) sum of squared deviations from center of n values,
 * continued from *ss and its compensation *kc */
static inline void column_kahan(const ST_double *RESTRICT x, const ST_double *RESTRICT w, ST_int n,
                                ST_double center, ST_double *ss_io, ST_double *kc_io)
{
    ST_double ss = *ss_io, kc = *kc_io;
    if (w)
    {
        for (ST_int i = 0; i < n; i++)
        {
            const ST_double d = x[i] - center;
            const ST_double sq = d * d;
            const ST_double wsq = w[i] * sq;
            const ST_double t = ss + wsq;
            kc += (ss - t) + wsq;
            ss = t;
        }
    }
    else
    {
        for (ST_int i = 0; i < n; i++)
        {
            const ST_double d = x[i] - center;
            const ST_double sq = d * d;
            const ST_double t = ss + sq;
            kc += (ss - t) + sq;
            ss = t;
        }
    }
    *ss_io = ss;
    *kc_io = kc;
}

/* Column k's compensated sum of squared deviations from center, in row order */
static ST_double column_ss(const creghdfe_columns *C, ST_int k, ST_double center)
{
    const ST_int N = C->N;
    ST_double buf[CREGHDFE_COLUMN_BLOCK];
    ST_double ss = 0.0, kc = 0.0;
    for (ST_int i0 = 0; i0 < N; i0 += CREGHDFE_COLUMN_BLOCK)
    {
        const ST_int n = N - i0 < CREGHDFE_COLUMN_BLOCK ? N - i0 : CREGHDFE_COLUMN_BLOCK;
        column_kahan(column_values(C, k, i0, n, buf), C->weights ? C->weights + i0 : NULL, n,
                     center, &ss, &kc);
    }
    return ss + kc;
}

/*
 * dst[0..n) = x * inv (x itself when !scale), and its squares added to the
 * eight chains s[] exactly as ctools_dot_unrolled(dst, dst) adds them: value
 * i to chain i mod 8, counting from the first row of the column, which the
 * caller guarantees by passing whole multiples of eight until the last
 * call. The values past the last multiple of eight are left to the caller
 * (*done receives how many were taken).
 */
static inline void column_store(const ST_double *RESTRICT x, ST_int n, ST_double inv, int scale,
                                ST_double *RESTRICT dst, ST_double *RESTRICT s, ST_int *done)
{
    const ST_int main = n - (n % 8);
    ST_double s0 = s[0], s1 = s[1], s2 = s[2], s3 = s[3];
    ST_double s4 = s[4], s5 = s[5], s6 = s[6], s7 = s[7];
    for (ST_int i = 0; i < main; i += 8)
    {
        const ST_double v0 = scale ? x[i] * inv : x[i];
        const ST_double v1 = scale ? x[i + 1] * inv : x[i + 1];
        const ST_double v2 = scale ? x[i + 2] * inv : x[i + 2];
        const ST_double v3 = scale ? x[i + 3] * inv : x[i + 3];
        const ST_double v4 = scale ? x[i + 4] * inv : x[i + 4];
        const ST_double v5 = scale ? x[i + 5] * inv : x[i + 5];
        const ST_double v6 = scale ? x[i + 6] * inv : x[i + 6];
        const ST_double v7 = scale ? x[i + 7] * inv : x[i + 7];
        dst[i] = v0;
        dst[i + 1] = v1;
        dst[i + 2] = v2;
        dst[i + 3] = v3;
        dst[i + 4] = v4;
        dst[i + 5] = v5;
        dst[i + 6] = v6;
        dst[i + 7] = v7;
        s0 += v0 * v0;
        s1 += v1 * v1;
        s2 += v2 * v2;
        s3 += v3 * v3;
        s4 += v4 * v4;
        s5 += v5 * v5;
        s6 += v6 * v6;
        s7 += v7 * v7;
    }
    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
    s[4] = s4;
    s[5] = s5;
    s[6] = s6;
    s[7] = s7;
    *done = main;
}

#ifdef CREGHDFE_STREAM_STORES
/* column_store with stores that do not read the destination first. odd: dst
 * is eight bytes off a 16-byte boundary, and the pairs are stored one value
 * later so that each still fills an aligned slot. The products and sums are
 * those of column_store, element by element. */
static inline void column_store_stream(const ST_double *RESTRICT x, ST_int n, ST_double inv,
                                       int scale, int odd, ST_double *RESTRICT dst,
                                       ST_double *RESTRICT s, ST_int *done)
{
    const ST_int main = n - (n % 8);
    const __m128d vinv = _mm_set1_pd(inv);
    __m128d a01 = _mm_loadu_pd(s), a23 = _mm_loadu_pd(s + 2);
    __m128d a45 = _mm_loadu_pd(s + 4), a67 = _mm_loadu_pd(s + 6);
    for (ST_int i = 0; i < main; i += 8)
    {
        __m128d v01 = _mm_loadu_pd(x + i), v23 = _mm_loadu_pd(x + i + 2);
        __m128d v45 = _mm_loadu_pd(x + i + 4), v67 = _mm_loadu_pd(x + i + 6);
        if (scale)
        {
            v01 = _mm_mul_pd(v01, vinv);
            v23 = _mm_mul_pd(v23, vinv);
            v45 = _mm_mul_pd(v45, vinv);
            v67 = _mm_mul_pd(v67, vinv);
        }
        if (odd)
        {
            _mm_stream_si64((long long *)(void *)(dst + i), _mm_cvtsi128_si64(_mm_castpd_si128(v01)));
            _mm_stream_pd(dst + i + 1, _mm_shuffle_pd(v01, v23, 1));
            _mm_stream_pd(dst + i + 3, _mm_shuffle_pd(v23, v45, 1));
            _mm_stream_pd(dst + i + 5, _mm_shuffle_pd(v45, v67, 1));
            _mm_stream_si64((long long *)(void *)(dst + i + 7),
                            _mm_cvtsi128_si64(_mm_castpd_si128(_mm_unpackhi_pd(v67, v67))));
        }
        else
        {
            _mm_stream_pd(dst + i, v01);
            _mm_stream_pd(dst + i + 2, v23);
            _mm_stream_pd(dst + i + 4, v45);
            _mm_stream_pd(dst + i + 6, v67);
        }
        a01 = _mm_add_pd(a01, _mm_mul_pd(v01, v01));
        a23 = _mm_add_pd(a23, _mm_mul_pd(v23, v23));
        a45 = _mm_add_pd(a45, _mm_mul_pd(v45, v45));
        a67 = _mm_add_pd(a67, _mm_mul_pd(v67, v67));
    }
    _mm_storeu_pd(s, a01);
    _mm_storeu_pd(s + 2, a23);
    _mm_storeu_pd(s + 4, a45);
    _mm_storeu_pd(s + 6, a67);
    *done = main;
}
#endif

/* Column k written to dst (and unscaled to raw), with the sum of squares of
 * what was written: the value ctools_dot_unrolled(dst, dst, N) returns */
static ST_double column_write(const creghdfe_columns *C, ST_int k, ST_double inv, int scale,
                              int stream, ST_double *RESTRICT dst, ST_double *RESTRICT raw)
{
    const ST_int N = C->N;
    ST_double buf[CREGHDFE_COLUMN_BLOCK];
    ST_double s[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef CREGHDFE_STREAM_STORES
    /* a destination that is not even 8-byte aligned is stored the ordinary way */
    if (((uintptr_t)dst & 7) != 0)
        stream = 0;
    const int odd = ((uintptr_t)dst & 15) != 0;
#else
    (void)stream;
#endif
    for (ST_int i0 = 0; i0 < N; i0 += CREGHDFE_COLUMN_BLOCK)
    {
        const ST_int n = N - i0 < CREGHDFE_COLUMN_BLOCK ? N - i0 : CREGHDFE_COLUMN_BLOCK;
        const ST_double *RESTRICT x = column_values(C, k, i0, n, buf);
        ST_int done = 0;
        if (raw)
            memcpy(raw + i0, x, (size_t)n * sizeof(ST_double));
#ifdef CREGHDFE_STREAM_STORES
        if (stream)
            column_store_stream(x, n, inv, scale, odd, dst + i0, s, &done);
        else
#endif
            column_store(x, n, inv, scale, dst + i0, s, &done);
        /* the last rows of the column, fewer than eight, join the first chain */
        for (ST_int i = done; i < n; i++)
        {
            const ST_double v = scale ? x[i] * inv : x[i];
            dst[i0 + i] = v;
            s[0] += v * v;
        }
    }
#ifdef CREGHDFE_STREAM_STORES
    if (stream)
        _mm_sfence();
#endif
    return (((s[0] + s[4]) + (s[1] + s[5])) + ((s[2] + s[6]) + (s[3] + s[7])));
}

void creghdfe_columns_build(const creghdfe_columns *C, const creghdfe_columns_opts *opts,
                            ST_double *data, ST_double *raw, ST_double *tss, ST_double *stdevs,
                            ST_double *norm2)
{
    const ST_int N = C->N, K = C->K;
    const ST_double *weights = C->weights;
    const ST_double df = (weights && opts->weight_type == 2) ? opts->sum_weights - 1.0
                                                             : (ST_double)(N - 1);
    const int stream = (size_t)N * (size_t)K * sizeof(ST_double) > CREGHDFE_STREAM_BYTES;
    ST_int k;
#pragma omp parallel for schedule(static) if (K > 1 && (size_t)N * (size_t)K > 200000)
    for (k = 0; k < K; k++)
    {
        ST_double *RESTRICT dst = data + (size_t)k * N;
        ST_double *RESTRICT raw_k = raw ? raw + (size_t)k * N : NULL;
        const ST_double mean = opts->means[k];
        ST_double ss;
        if (opts->use_quad)
        {
            /* the double-double sums take the column in one piece: copy it,
             * then scale the copy where it lies */
            ST_double buf[CREGHDFE_COLUMN_BLOCK];
            for (ST_int i0 = 0; i0 < N; i0 += CREGHDFE_COLUMN_BLOCK)
            {
                const ST_int n = N - i0 < CREGHDFE_COLUMN_BLOCK ? N - i0 : CREGHDFE_COLUMN_BLOCK;
                memcpy(dst + i0, column_values(C, k, i0, n, buf), (size_t)n * sizeof(ST_double));
            }
            if (raw_k)
                memcpy(raw_k, dst, (size_t)N * sizeof(ST_double));
            ss = weights ? dd_sum_sq_dev_weighted(dst, mean, weights, N)
                         : dd_sum_sq_dev(dst, mean, N);
        }
        else
            ss = column_ss(C, k, mean);
        tss[k] = ss;
        stdevs[k] = sqrt(ss / df);
        if (stdevs[k] < 1e-30)
            stdevs[k] = 1.0;
        /* Without a constant the total sum of squares is uncentered, as in
         * reghdfe (the standard deviations stay centered) */
        if (opts->uncentered)
            tss[k] = column_ss(C, k, 0.0);
        /* An arbitrary response multiplier rounds small identifying
         * components beside a large fitted component before OLS ever sees
         * them. A power of two preserves the response's bits (away from
         * under/overflow) and keeps its standard deviation in [1,2).
         * Regressor scaling and rank decisions are unchanged. */
        if (k == 0 && opts->standardize && isfinite(stdevs[0]) && stdevs[0] > 0.0)
        {
            int exponent;
            frexp(stdevs[0], &exponent);
            stdevs[0] = ldexp(1.0, exponent - 1);
        }
        const ST_double inv = 1.0 / stdevs[k];
        if (opts->use_quad)
        {
            if (opts->standardize)
            {
#pragma omp simd
                for (ST_int i = 0; i < N; i++)
                    dst[i] *= inv;
            }
            norm2[k] = fast_dot(dst, dst, N);
        }
        else
            norm2[k] = column_write(C, k, inv, opts->standardize != 0, stream, dst, raw_k);
    }
}

/* ========================================================================
 * Collinearity and OLS on the partialled columns
 * ======================================================================== */

ST_retcode creghdfe_collinear(ST_double *data, ST_int N, ST_int K, const ST_double *norm2,
                              const ST_double *count_w, const ST_double *weights, ST_int weight_type,
                              ST_double collinear_tol, ST_int *fe_flagged,
                              ST_int *is_collinear, ST_int *depvar_explained, ST_double *xtx,
                              ST_double *xty_all, ST_double *scratch)
{
    const ST_int K_x = K - 1;
    ST_int k, j;
    if (K_x > 0)
    {
        memset(is_collinear, 0, (size_t)K_x * sizeof(ST_int));
        memset(fe_flagged, 0, (size_t)K_x * sizeof(ST_int));
    }

    /* Columns collinear with the fixed effects, by reghdfe's rule
     * (check_collinear_with_fe): the partialled sum of squares is at most
     * collinear_tol = min(1e-6, tolerance()/10) times the uncentered sum of
     * squares before partialling (both unweighted, so the rule does not
     * depend on the units or on standardization). Such a regressor is
     * omitted. A dependent variable that passes the rule is perfectly
     * explained by the effects; it is zeroed, as reghdfe does, so its
     * coefficients, rss and tss_within are exactly zero. */
    *depvar_explained = 0;
    ST_int num_fe_collinear = 0;
    /* The partialled sums of squares: one independent sum per column, each in
     * row order, so the columns run in parallel */
    ST_double *ss_all = (ST_double *)malloc((size_t)K * sizeof(ST_double));
    if (!ss_all)
        return 920;
#pragma omp parallel for schedule(static) if (K > 1 && (size_t)N * (size_t)K > 200000)
    for (k = 0; k < K; k++)
    {
        const ST_double *col = &data[(size_t)k * N];
        ST_double ss;
        if (count_w)
        {
            ss = 0.0;
            for (ST_int i = 0; i < N; i++)
                ss += count_w[i] * col[i] * col[i];
        }
        else
            ss = fast_dot(col, col, N);
        ss_all[k] = ss;
    }
    for (k = 0; k < K; k++)
    {
        ST_double *col = &data[(size_t)k * N];
        ST_double ss = ss_all[k];
        ST_double ratio = ss / norm2[k];
        ST_int flagged = (ratio <= collinear_tol);
        if (k == 0)
        {
            if (flagged)
            {
                memset(col, 0, (size_t)N * sizeof(ST_double));
                *depvar_explained = 1;
            }
            continue;
        }
        if (flagged)
        {
            is_collinear[k - 1] = 1;
            num_fe_collinear++;
        }
        fe_flagged[k - 1] = flagged;
    }
    free(ss_all);

    /* The remaining collinearity as reghdfe_rmcoll: Mata's invsym() of the
     * (weighted) cross products of the partialled regressors, sweeping them
     * in order, and invsym() pivoting on the largest diagonal instead when
     * that drops a different number of them. A regressor is collinear when
     * its diagonal of the generalized inverse is zero. */
    if (K_x > 0)
    {
        if (weights)
            compute_xtx_xty_weighted(data, weights, weight_type, N, K, xtx, xty_all);
        else
            compute_xtx_xty(data, N, K, xtx, xty_all);
        ST_double *xx_rm = xtx;
        if (num_fe_collinear > 0)
        {
            /* Those columns leave the cross products, as reghdfe trims them */
            xx_rm = scratch;
            memcpy(xx_rm, xtx, (size_t)K_x * K_x * sizeof(ST_double));
            for (k = 0; k < K_x; k++)
            {
                if (!is_collinear[k])
                    continue;
                for (j = 0; j < K_x; j++)
                {
                    xx_rm[(size_t)k * K_x + j] = 0.0;
                    xx_rm[(size_t)j * K_x + k] = 0.0;
                }
            }
        }
        ST_double *ginv = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
        if (!ginv)
            return 920;
        ST_int rank_ordered = 0, rank_pivoted = 0;
        ST_int inv_rc = ctools_invsym_ordered(xx_rm, K_x, ginv, &rank_ordered);
        if (!inv_rc && rank_ordered < K_x - num_fe_collinear)
        {
            ST_double *ginv_alt = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
            if (!ginv_alt)
            {
                free(ginv);
                return 920;
            }
            inv_rc = ctools_invsym(xx_rm, K_x, ginv_alt, &rank_pivoted);
            if (!inv_rc && rank_pivoted != rank_ordered)
            {
                free(ginv);
                ginv = ginv_alt;
            }
            else
            {
                free(ginv_alt);
            }
        }
        if (inv_rc)
        {
            free(ginv);
            return 920;
        }
        for (k = 0; k < K_x; k++)
        {
            if (ginv[(size_t)k * K_x + k] == 0.0)
                is_collinear[k] = 1;
        }
        free(ginv);
    }
    return 0;
}

ST_retcode creghdfe_ols(const creghdfe_ols_in *in, creghdfe_ols_out *out)
{
    const ST_int N = in->N, K_x = in->K - 1;
    const ST_double *weights = in->weights;
    ST_int i, j, k, idx;
    memset(out, 0, sizeof(*out));

    ST_int num_collinear = 0;
    for (k = 0; k < K_x; k++)
    {
        if (in->is_collinear[k])
            num_collinear++;
    }
    const ST_int K_keep = K_x - num_collinear;
    out->K_keep = K_keep;

    /* K_keep == 0 (every regressor omitted) continues through the same path,
     * as reghdfe does: the model is then the fixed effects and the constant.
     * K_alloc keeps the K_keep-sized allocations nonempty in that case. */
    const size_t K_alloc = (K_keep > 0) ? (size_t)K_keep : 1;

    /* Build index of non-collinear variables */
    out->keep_idx = (ST_int *)malloc(K_alloc * sizeof(ST_int));
    if (!out->keep_idx)
        return 920;
    idx = 0;
    for (k = 0; k < K_x; k++)
    {
        if (!in->is_collinear[k])
            out->keep_idx[idx++] = k + 1; /* +1 to skip y */
    }

    const ST_int K_with_cons = K_keep + 1;
    /* y and the kept regressors. When nothing is omitted this is the data
     * matrix itself and no copy is made. */
    const int keep_is_data = (K_keep == K_x);
    out->data_keep = keep_is_data
                         ? in->data
                         : (ST_double *)malloc((size_t)N * (K_keep + 1) * sizeof(ST_double));
    out->xtx_keep = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    out->xty_keep = (ST_double *)malloc(K_alloc * sizeof(ST_double));
    out->beta_keep = (ST_double *)malloc(K_with_cons * sizeof(ST_double));
    out->inv_xx_keep = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    out->means_x = (ST_double *)malloc(K_alloc * sizeof(ST_double));
    if (!out->data_keep || !out->xtx_keep || !out->xty_keep || !out->beta_keep ||
        !out->inv_xx_keep || !out->means_x)
        return 920;
    ST_double *data_keep = out->data_keep, *xtx_keep = out->xtx_keep;
    ST_double *xty_keep = out->xty_keep, *beta_keep = out->beta_keep;
    ST_double *means_x = out->means_x;
    const ST_int *keep_idx = out->keep_idx;

    /* Copy y and non-collinear X */
    if (!keep_is_data)
    {
        memcpy(data_keep, in->data, (size_t)N * sizeof(ST_double));
        for (k = 0; k < K_keep; k++)
        {
            memcpy(data_keep + (size_t)(k + 1) * N, in->data + (size_t)keep_idx[k] * N,
                   (size_t)N * sizeof(ST_double));
        }
    }
    for (k = 0; k < K_keep; k++)
    {
        if (in->standardize)
            means_x[k] = in->means[keep_idx[k]] / in->stdevs[keep_idx[k]];
        else
            means_x[k] = in->means[keep_idx[k]];
    }

    /* Compute X'X and X'y on non-collinear data (weighted if using weights);
     * when nothing is omitted they are the cross products computed before */
    ST_double tss_within;
    const int reuse_xx = keep_is_data && K_keep > 0;
    if (reuse_xx)
    {
        memcpy(xtx_keep, in->xtx, (size_t)K_keep * K_keep * sizeof(ST_double));
        memcpy(xty_keep, in->xty_all, (size_t)K_keep * sizeof(ST_double));
    }
    if (weights)
    {
        if (!reuse_xx)
            compute_xtx_xty_weighted(data_keep, weights, in->weight_type, N, K_keep + 1, xtx_keep,
                                     xty_keep);
        /* Weighted TSS_within = sum(w_i * y_i^2)
         * Use quad-precision or Kahan summation based on option. */
        if (in->use_quad)
        {
            tss_within = dd_sum_sq_weighted(data_keep, weights, N);
        }
        else
        {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++)
            {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += weights[idx] * sq;
            }
        }
    }
    else
    {
        if (!reuse_xx)
            compute_xtx_xty(data_keep, N, K_keep + 1, xtx_keep, xty_keep);
        /* Compute tss_within: use quad-precision or Kahan summation. */
        if (in->use_quad)
        {
            tss_within = dd_sum_sq(data_keep, N);
        }
        else
        {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++)
            {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += sq;
            }
        }
    }
    out->tss_within = tss_within;

    /* Cholesky and invert */
    ST_double *inv_xx_x = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    if (!inv_xx_x)
        return 920;
    if (K_keep > 0)
    {
        memcpy(inv_xx_x, xtx_keep, (size_t)K_keep * K_keep * sizeof(ST_double));
        if (ctools_cholesky(inv_xx_x, K_keep) != 0)
        {
            free(inv_xx_x);
            return 198;
        }
        stata_retcode inverse_status_0 = ctools_invert_from_cholesky(inv_xx_x, K_keep, inv_xx_x);
        if (inverse_status_0 != STATA_OK)
        {
            free(inv_xx_x);
            return ctools_stata_rc(inverse_status_0);
        }
    }

    /* Compute beta */
    for (i = 0; i < K_keep; i++)
    {
        beta_keep[i] = 0.0;
        for (j = 0; j < K_keep; j++)
            beta_keep[i] += inv_xx_x[i * K_keep + j] * xty_keep[j];
    }

    /* Compute RSS */
    ST_double rss = tss_within;
    for (k = 0; k < K_keep; k++)
        rss -= beta_keep[k] * xty_keep[k];
    out->rss = rss;

    /* Extend inv_xx with the constant (block partition formula). For
     * fweight the 1'W1 term is sum(weights). */
    ST_double N_corner = (in->weight_type == 2 && weights) ? in->sum_weights : (ST_double)N;
    ctools_inverse_add_constant(inv_xx_x, K_keep, means_x, N_corner, out->inv_xx_keep);

    /* Compute constant */
    ST_double y_mean = in->means[0];
    ST_double xb_mean = 0.0;
    for (k = 0; k < K_keep; k++)
        xb_mean += means_x[k] * beta_keep[k];
    beta_keep[K_keep] = y_mean - xb_mean;

    free(inv_xx_x);
    return 0;
}

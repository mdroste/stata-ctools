/*
 * ctools_simd.h
 *
 * SIMD (AVX2/NEON) abstractions for ctools.
 * Provides portable SIMD utilities for vectorized operations on double arrays.
 *
 * Supported architectures:
 *   - AVX2 (x86_64): 4 doubles per operation (256-bit)
 *   - NEON (ARM64): 2 doubles per operation (128-bit)
 *
 * Usage:
 *   #include "ctools_simd.h"
 *
 *   #if CTOOLS_HAS_AVX2
 *       // AVX2 code path
 *   #elif CTOOLS_HAS_NEON
 *       // NEON code path
 *   #endif
 *   // scalar fallback
 *
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_SIMD_H
#define CTOOLS_SIMD_H

#include "stplugin.h"
#include <float.h>

/* ============================================================================
 * Platform Detection and Include
 * ============================================================================ */

#if defined(__AVX2__)
    #include <immintrin.h>
    #define CTOOLS_HAS_AVX2 1
    #define CTOOLS_HAS_SSE2 1       /* AVX2 implies SSE2 */
    #define CTOOLS_HAS_NEON 0
    #define CTOOLS_SIMD_WIDTH 4     /* doubles per vector */
    #define CTOOLS_SIMD_BYTES 32    /* bytes per vector */
#elif defined(__SSE2__)
    #include <emmintrin.h>
    #define CTOOLS_HAS_AVX2 0
    #define CTOOLS_HAS_SSE2 1
    #define CTOOLS_HAS_NEON 0
    #define CTOOLS_SIMD_WIDTH 2     /* doubles per vector */
    #define CTOOLS_SIMD_BYTES 16    /* bytes per vector */
#elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
    #define CTOOLS_HAS_AVX2 0
    #define CTOOLS_HAS_SSE2 0
    #define CTOOLS_HAS_NEON 1
    #define CTOOLS_SIMD_WIDTH 2     /* doubles per vector */
    #define CTOOLS_SIMD_BYTES 16    /* bytes per vector */
#else
    #define CTOOLS_HAS_AVX2 0
    #define CTOOLS_HAS_SSE2 0
    #define CTOOLS_HAS_NEON 0
    #define CTOOLS_SIMD_WIDTH 1
    #define CTOOLS_SIMD_BYTES 8
#endif

/* Minimum elements to use SIMD path (below this, scalar is faster) */
#define CTOOLS_SIMD_THRESHOLD 8

/* ============================================================================
 * AVX2 Horizontal Reduction Utilities
 * ============================================================================ */

#if CTOOLS_HAS_AVX2

/*
 * ctools_hsum_avx2 - Horizontal sum of 4 doubles in AVX2 vector
 *
 * Uses tree reduction for numerical stability:
 *   result = (v[0] + v[2]) + (v[1] + v[3])
 */
static inline double ctools_hsum_avx2(__m256d v)
{
    /* Extract high 128 bits and add to low 128 bits */
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    __m128d vsum  = _mm_add_pd(vlow, vhigh);

    /* Horizontal add within 128-bit register */
    __m128d vshuf = _mm_shuffle_pd(vsum, vsum, 1);
    __m128d vresult = _mm_add_sd(vsum, vshuf);

    return _mm_cvtsd_f64(vresult);
}

/*
 * ctools_hmin_avx2 - Horizontal minimum of 4 doubles
 */
static inline double ctools_hmin_avx2(__m256d v)
{
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    __m128d vmin  = _mm_min_pd(vlow, vhigh);
    __m128d vshuf = _mm_shuffle_pd(vmin, vmin, 1);
    __m128d vresult = _mm_min_sd(vmin, vshuf);
    return _mm_cvtsd_f64(vresult);
}

/*
 * ctools_hmax_avx2 - Horizontal maximum of 4 doubles
 */
static inline double ctools_hmax_avx2(__m256d v)
{
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    __m128d vmax  = _mm_max_pd(vlow, vhigh);
    __m128d vshuf = _mm_shuffle_pd(vmax, vmax, 1);
    __m128d vresult = _mm_max_sd(vmax, vshuf);
    return _mm_cvtsd_f64(vresult);
}

/*
 * ctools_valid_mask_avx2 - Create mask for non-missing values
 *
 * Returns mask where lane is all-1s if val < SV_missval, all-0s otherwise.
 * Use with _mm256_blendv_pd() for conditional operations.
 */
static inline __m256d ctools_valid_mask_avx2(__m256d v, __m256d vmiss)
{
    return _mm256_cmp_pd(v, vmiss, _CMP_LT_OQ);
}

/*
 * ctools_count_valid_avx2 - Count non-missing values in vector
 *
 * Returns count of lanes where val < SV_missval (0-4).
 */
static inline int ctools_count_valid_avx2(__m256d mask)
{
    /* movemask returns 4-bit mask, popcount gives count */
    int bits = _mm256_movemask_pd(mask);
    return __builtin_popcount(bits);
}

#endif /* CTOOLS_HAS_AVX2 */

/* ============================================================================
 * NEON Horizontal Reduction Utilities
 * ============================================================================ */

#if CTOOLS_HAS_NEON

/*
 * ctools_hsum_neon - Horizontal sum of 2 doubles in NEON vector
 */
static inline double ctools_hsum_neon(float64x2_t v)
{
    return vgetq_lane_f64(v, 0) + vgetq_lane_f64(v, 1);
}

/*
 * ctools_hmin_neon - Horizontal minimum of 2 doubles
 */
static inline double ctools_hmin_neon(float64x2_t v)
{
    double a = vgetq_lane_f64(v, 0);
    double b = vgetq_lane_f64(v, 1);
    return (a < b) ? a : b;
}

/*
 * ctools_hmax_neon - Horizontal maximum of 2 doubles
 */
static inline double ctools_hmax_neon(float64x2_t v)
{
    double a = vgetq_lane_f64(v, 0);
    double b = vgetq_lane_f64(v, 1);
    return (a > b) ? a : b;
}

/*
 * ctools_valid_mask_neon - Create mask for non-missing values
 *
 * Returns mask where lane is all-1s if val < SV_missval, all-0s otherwise.
 */
static inline uint64x2_t ctools_valid_mask_neon(float64x2_t v, float64x2_t vmiss)
{
    return vcltq_f64(v, vmiss);
}

/*
 * ctools_count_valid_neon - Count non-missing values in vector
 *
 * Returns count of lanes where val < SV_missval (0-2).
 */
static inline int ctools_count_valid_neon(uint64x2_t mask)
{
    uint64_t bits[2];
    vst1q_u64(bits, mask);
    return (bits[0] ? 1 : 0) + (bits[1] ? 1 : 0);
}

#endif /* CTOOLS_HAS_NEON */

/* ============================================================================
 * Portable SIMD Operations
 * ============================================================================ */

/*
 * ctools_simd_clamp - Clamp values to [lower, upper] with SIMD
 *
 * Processes array in place, skipping missing values.
 * Uses AVX2/NEON when available, scalar fallback otherwise.
 */
static inline void ctools_simd_clamp(double *data, size_t n,
                                     double lower, double upper,
                                     double miss)
{
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d vlower = _mm256_set1_pd(lower);
        __m256d vupper = _mm256_set1_pd(upper);
        __m256d vmiss  = _mm256_set1_pd(miss);

        for (; i + 4 <= n; i += 4) {
            __m256d v = _mm256_loadu_pd(&data[i]);
            __m256d valid = _mm256_cmp_pd(v, vmiss, _CMP_LT_OQ);

            /* Clamp: max(lower, min(upper, v)) */
            __m256d clamped = _mm256_max_pd(vlower, _mm256_min_pd(vupper, v));

            /* Blend: use clamped where valid, original where missing */
            __m256d result = _mm256_blendv_pd(v, clamped, valid);
            _mm256_storeu_pd(&data[i], result);
        }
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t vlower = vdupq_n_f64(lower);
        float64x2_t vupper = vdupq_n_f64(upper);
        float64x2_t vmiss  = vdupq_n_f64(miss);

        for (; i + 2 <= n; i += 2) {
            float64x2_t v = vld1q_f64(&data[i]);
            uint64x2_t valid = vcltq_f64(v, vmiss);

            /* Clamp: max(lower, min(upper, v)) */
            float64x2_t clamped = vmaxq_f64(vlower, vminq_f64(vupper, v));

            /* Blend: use clamped where valid, original where missing */
            float64x2_t result = vbslq_f64(valid, clamped, v);
            vst1q_f64(&data[i], result);
        }
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        double val = data[i];
        if (val < miss) {
            if (val < lower) data[i] = lower;
            else if (val > upper) data[i] = upper;
        }
    }
}

/*
 * ctools_simd_replace_oob - Replace out-of-bounds values with missing
 *
 * Processes array in place, replacing values outside [lower, upper] with miss.
 * Used for trimming operations.
 */
static inline void ctools_simd_replace_oob(double *data, size_t n,
                                           double lower, double upper,
                                           double miss)
{
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d vlower = _mm256_set1_pd(lower);
        __m256d vupper = _mm256_set1_pd(upper);
        __m256d vmiss  = _mm256_set1_pd(miss);

        for (; i + 4 <= n; i += 4) {
            __m256d v = _mm256_loadu_pd(&data[i]);

            /* Check: valid AND in_bounds */
            __m256d valid = _mm256_cmp_pd(v, vmiss, _CMP_LT_OQ);
            __m256d ge_lower = _mm256_cmp_pd(v, vlower, _CMP_GE_OQ);
            __m256d le_upper = _mm256_cmp_pd(v, vupper, _CMP_LE_OQ);
            __m256d in_bounds = _mm256_and_pd(ge_lower, le_upper);
            __m256d keep = _mm256_and_pd(valid, in_bounds);

            /* Replace out-of-bounds with missing */
            __m256d result = _mm256_blendv_pd(vmiss, v, keep);
            _mm256_storeu_pd(&data[i], result);
        }
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t vlower = vdupq_n_f64(lower);
        float64x2_t vupper = vdupq_n_f64(upper);
        float64x2_t vmiss  = vdupq_n_f64(miss);

        for (; i + 2 <= n; i += 2) {
            float64x2_t v = vld1q_f64(&data[i]);

            /* Check: valid AND in_bounds */
            uint64x2_t valid = vcltq_f64(v, vmiss);
            uint64x2_t ge_lower = vcgeq_f64(v, vlower);
            uint64x2_t le_upper = vcleq_f64(v, vupper);
            uint64x2_t in_bounds = vandq_u64(ge_lower, le_upper);
            uint64x2_t keep = vandq_u64(valid, in_bounds);

            /* Replace out-of-bounds with missing */
            float64x2_t result = vbslq_f64(keep, v, vmiss);
            vst1q_f64(&data[i], result);
        }
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        double val = data[i];
        if (val < miss) {
            if (val < lower || val > upper) data[i] = miss;
        }
    }
}

/*
 * ctools_simd_dot - SIMD-accelerated dot product
 *
 * Computes sum(x[i] * y[i]) using AVX2 FMA when available.
 */
static inline double ctools_simd_dot(const double *x, const double *y, size_t n)
{
    double sum = 0.0;
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d vsum0 = _mm256_setzero_pd();
        __m256d vsum1 = _mm256_setzero_pd();

        /* 8-way unroll with 2 accumulators */
        for (; i + 8 <= n; i += 8) {
            __m256d vx0 = _mm256_loadu_pd(&x[i]);
            __m256d vy0 = _mm256_loadu_pd(&y[i]);
            __m256d vx1 = _mm256_loadu_pd(&x[i + 4]);
            __m256d vy1 = _mm256_loadu_pd(&y[i + 4]);

            #ifdef __FMA__
            vsum0 = _mm256_fmadd_pd(vx0, vy0, vsum0);
            vsum1 = _mm256_fmadd_pd(vx1, vy1, vsum1);
            #else
            vsum0 = _mm256_add_pd(vsum0, _mm256_mul_pd(vx0, vy0));
            vsum1 = _mm256_add_pd(vsum1, _mm256_mul_pd(vx1, vy1));
            #endif
        }

        /* 4-element cleanup */
        for (; i + 4 <= n; i += 4) {
            __m256d vx = _mm256_loadu_pd(&x[i]);
            __m256d vy = _mm256_loadu_pd(&y[i]);
            #ifdef __FMA__
            vsum0 = _mm256_fmadd_pd(vx, vy, vsum0);
            #else
            vsum0 = _mm256_add_pd(vsum0, _mm256_mul_pd(vx, vy));
            #endif
        }

        __m256d vtotal = _mm256_add_pd(vsum0, vsum1);
        sum = ctools_hsum_avx2(vtotal);
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t vsum0 = vdupq_n_f64(0.0);
        float64x2_t vsum1 = vdupq_n_f64(0.0);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vx0 = vld1q_f64(&x[i]);
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vx1 = vld1q_f64(&x[i + 2]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);

            vsum0 = vfmaq_f64(vsum0, vx0, vy0);
            vsum1 = vfmaq_f64(vsum1, vx1, vy1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vx = vld1q_f64(&x[i]);
            float64x2_t vy = vld1q_f64(&y[i]);
            vsum0 = vfmaq_f64(vsum0, vx, vy);
        }

        float64x2_t vtotal = vaddq_f64(vsum0, vsum1);
        sum = ctools_hsum_neon(vtotal);
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        sum += x[i] * y[i];
    }

    return sum;
}

/*
 * ctools_simd_dot_weighted - SIMD-accelerated weighted dot product
 *
 * Computes sum(w[i] * x[i] * y[i]) using AVX2 FMA when available.
 */
static inline double ctools_simd_dot_weighted(const double *w,
                                              const double *x,
                                              const double *y,
                                              size_t n)
{
    double sum = 0.0;
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d vsum0 = _mm256_setzero_pd();
        __m256d vsum1 = _mm256_setzero_pd();

        for (; i + 8 <= n; i += 8) {
            __m256d vw0 = _mm256_loadu_pd(&w[i]);
            __m256d vx0 = _mm256_loadu_pd(&x[i]);
            __m256d vy0 = _mm256_loadu_pd(&y[i]);
            __m256d vw1 = _mm256_loadu_pd(&w[i + 4]);
            __m256d vx1 = _mm256_loadu_pd(&x[i + 4]);
            __m256d vy1 = _mm256_loadu_pd(&y[i + 4]);

            __m256d vwx0 = _mm256_mul_pd(vw0, vx0);
            __m256d vwx1 = _mm256_mul_pd(vw1, vx1);

            #ifdef __FMA__
            vsum0 = _mm256_fmadd_pd(vwx0, vy0, vsum0);
            vsum1 = _mm256_fmadd_pd(vwx1, vy1, vsum1);
            #else
            vsum0 = _mm256_add_pd(vsum0, _mm256_mul_pd(vwx0, vy0));
            vsum1 = _mm256_add_pd(vsum1, _mm256_mul_pd(vwx1, vy1));
            #endif
        }

        for (; i + 4 <= n; i += 4) {
            __m256d vw = _mm256_loadu_pd(&w[i]);
            __m256d vx = _mm256_loadu_pd(&x[i]);
            __m256d vy = _mm256_loadu_pd(&y[i]);
            __m256d vwx = _mm256_mul_pd(vw, vx);

            #ifdef __FMA__
            vsum0 = _mm256_fmadd_pd(vwx, vy, vsum0);
            #else
            vsum0 = _mm256_add_pd(vsum0, _mm256_mul_pd(vwx, vy));
            #endif
        }

        __m256d vtotal = _mm256_add_pd(vsum0, vsum1);
        sum = ctools_hsum_avx2(vtotal);
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t vsum0 = vdupq_n_f64(0.0);
        float64x2_t vsum1 = vdupq_n_f64(0.0);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vw0 = vld1q_f64(&w[i]);
            float64x2_t vx0 = vld1q_f64(&x[i]);
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vw1 = vld1q_f64(&w[i + 2]);
            float64x2_t vx1 = vld1q_f64(&x[i + 2]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);

            float64x2_t vwx0 = vmulq_f64(vw0, vx0);
            float64x2_t vwx1 = vmulq_f64(vw1, vx1);

            vsum0 = vfmaq_f64(vsum0, vwx0, vy0);
            vsum1 = vfmaq_f64(vsum1, vwx1, vy1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vw = vld1q_f64(&w[i]);
            float64x2_t vx = vld1q_f64(&x[i]);
            float64x2_t vy = vld1q_f64(&y[i]);
            float64x2_t vwx = vmulq_f64(vw, vx);
            vsum0 = vfmaq_f64(vsum0, vwx, vy);
        }

        float64x2_t vtotal = vaddq_f64(vsum0, vsum1);
        sum = ctools_hsum_neon(vtotal);
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        sum += w[i] * x[i] * y[i];
    }

    return sum;
}

/*
 * ctools_simd_axpy - SIMD-accelerated y += alpha * x
 *
 * Performs y[i] += alpha * x[i] for all i.
 */
static inline void ctools_simd_axpy(double *y, const double *x,
                                    double alpha, size_t n)
{
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d valpha = _mm256_set1_pd(alpha);

        for (; i + 8 <= n; i += 8) {
            __m256d vy0 = _mm256_loadu_pd(&y[i]);
            __m256d vx0 = _mm256_loadu_pd(&x[i]);
            __m256d vy1 = _mm256_loadu_pd(&y[i + 4]);
            __m256d vx1 = _mm256_loadu_pd(&x[i + 4]);

            #ifdef __FMA__
            vy0 = _mm256_fmadd_pd(valpha, vx0, vy0);
            vy1 = _mm256_fmadd_pd(valpha, vx1, vy1);
            #else
            vy0 = _mm256_add_pd(vy0, _mm256_mul_pd(valpha, vx0));
            vy1 = _mm256_add_pd(vy1, _mm256_mul_pd(valpha, vx1));
            #endif

            _mm256_storeu_pd(&y[i], vy0);
            _mm256_storeu_pd(&y[i + 4], vy1);
        }

        for (; i + 4 <= n; i += 4) {
            __m256d vy = _mm256_loadu_pd(&y[i]);
            __m256d vx = _mm256_loadu_pd(&x[i]);

            #ifdef __FMA__
            vy = _mm256_fmadd_pd(valpha, vx, vy);
            #else
            vy = _mm256_add_pd(vy, _mm256_mul_pd(valpha, vx));
            #endif

            _mm256_storeu_pd(&y[i], vy);
        }
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t valpha = vdupq_n_f64(alpha);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vx0 = vld1q_f64(&x[i]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);
            float64x2_t vx1 = vld1q_f64(&x[i + 2]);

            vy0 = vfmaq_f64(vy0, valpha, vx0);
            vy1 = vfmaq_f64(vy1, valpha, vx1);

            vst1q_f64(&y[i], vy0);
            vst1q_f64(&y[i + 2], vy1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vy = vld1q_f64(&y[i]);
            float64x2_t vx = vld1q_f64(&x[i]);
            vy = vfmaq_f64(vy, valpha, vx);
            vst1q_f64(&y[i], vy);
        }
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

/*
 * ctools_simd_axpby - SIMD-accelerated y = alpha * x + beta * y
 *
 * Performs y[i] = alpha * x[i] + beta * y[i] for all i.
 */
static inline void ctools_simd_axpby(double *y, const double *x,
                                     double alpha, double beta, size_t n)
{
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d valpha = _mm256_set1_pd(alpha);
        __m256d vbeta  = _mm256_set1_pd(beta);

        for (; i + 4 <= n; i += 4) {
            __m256d vy = _mm256_loadu_pd(&y[i]);
            __m256d vx = _mm256_loadu_pd(&x[i]);

            __m256d vax = _mm256_mul_pd(valpha, vx);
            #ifdef __FMA__
            vy = _mm256_fmadd_pd(vbeta, vy, vax);
            #else
            vy = _mm256_add_pd(vax, _mm256_mul_pd(vbeta, vy));
            #endif

            _mm256_storeu_pd(&y[i], vy);
        }
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t valpha = vdupq_n_f64(alpha);
        float64x2_t vbeta  = vdupq_n_f64(beta);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vx0 = vld1q_f64(&x[i]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);
            float64x2_t vx1 = vld1q_f64(&x[i + 2]);

            float64x2_t vax0 = vmulq_f64(valpha, vx0);
            float64x2_t vax1 = vmulq_f64(valpha, vx1);
            vy0 = vfmaq_f64(vax0, vbeta, vy0);
            vy1 = vfmaq_f64(vax1, vbeta, vy1);

            vst1q_f64(&y[i], vy0);
            vst1q_f64(&y[i + 2], vy1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vy = vld1q_f64(&y[i]);
            float64x2_t vx = vld1q_f64(&x[i]);

            float64x2_t vax = vmulq_f64(valpha, vx);
            vy = vfmaq_f64(vax, vbeta, vy);

            vst1q_f64(&y[i], vy);
        }
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}

/*
 * ctools_simd_fused_axpy2_dot - Fused y -= alpha*u, r -= alpha*v, return r·r
 *
 * Performs three operations in a single pass over N doubles:
 *   y[i] -= alpha * u[i]
 *   r[i] -= alpha * v[i]
 *   sum  += r[i] * r[i]   (after update)
 *
 * Returns the dot product r·r (after the axpy update).
 * Saves 2 full passes over N compared to separate axpy+axpy+dot.
 */
static inline double ctools_simd_fused_axpy2_dot(double *y, const double *u,
                                                  double *r, const double *v,
                                                  double alpha, size_t n)
{
    double dot = 0.0;
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d valpha = _mm256_set1_pd(alpha);
        __m256d vdot0 = _mm256_setzero_pd();
        __m256d vdot1 = _mm256_setzero_pd();

        for (; i + 8 <= n; i += 8) {
            /* Load */
            __m256d vy0 = _mm256_loadu_pd(&y[i]);
            __m256d vu0 = _mm256_loadu_pd(&u[i]);
            __m256d vr0 = _mm256_loadu_pd(&r[i]);
            __m256d vv0 = _mm256_loadu_pd(&v[i]);
            __m256d vy1 = _mm256_loadu_pd(&y[i + 4]);
            __m256d vu1 = _mm256_loadu_pd(&u[i + 4]);
            __m256d vr1 = _mm256_loadu_pd(&r[i + 4]);
            __m256d vv1 = _mm256_loadu_pd(&v[i + 4]);

            /* y -= alpha * u */
            #ifdef __FMA__
            vy0 = _mm256_fnmadd_pd(valpha, vu0, vy0);
            vy1 = _mm256_fnmadd_pd(valpha, vu1, vy1);
            #else
            vy0 = _mm256_sub_pd(vy0, _mm256_mul_pd(valpha, vu0));
            vy1 = _mm256_sub_pd(vy1, _mm256_mul_pd(valpha, vu1));
            #endif

            /* r -= alpha * v */
            #ifdef __FMA__
            vr0 = _mm256_fnmadd_pd(valpha, vv0, vr0);
            vr1 = _mm256_fnmadd_pd(valpha, vv1, vr1);
            #else
            vr0 = _mm256_sub_pd(vr0, _mm256_mul_pd(valpha, vv0));
            vr1 = _mm256_sub_pd(vr1, _mm256_mul_pd(valpha, vv1));
            #endif

            /* Store */
            _mm256_storeu_pd(&y[i], vy0);
            _mm256_storeu_pd(&y[i + 4], vy1);
            _mm256_storeu_pd(&r[i], vr0);
            _mm256_storeu_pd(&r[i + 4], vr1);

            /* Accumulate r·r */
            #ifdef __FMA__
            vdot0 = _mm256_fmadd_pd(vr0, vr0, vdot0);
            vdot1 = _mm256_fmadd_pd(vr1, vr1, vdot1);
            #else
            vdot0 = _mm256_add_pd(vdot0, _mm256_mul_pd(vr0, vr0));
            vdot1 = _mm256_add_pd(vdot1, _mm256_mul_pd(vr1, vr1));
            #endif
        }

        for (; i + 4 <= n; i += 4) {
            __m256d vy = _mm256_loadu_pd(&y[i]);
            __m256d vu = _mm256_loadu_pd(&u[i]);
            __m256d vr = _mm256_loadu_pd(&r[i]);
            __m256d vv = _mm256_loadu_pd(&v[i]);

            #ifdef __FMA__
            vy = _mm256_fnmadd_pd(valpha, vu, vy);
            vr = _mm256_fnmadd_pd(valpha, vv, vr);
            #else
            vy = _mm256_sub_pd(vy, _mm256_mul_pd(valpha, vu));
            vr = _mm256_sub_pd(vr, _mm256_mul_pd(valpha, vv));
            #endif

            _mm256_storeu_pd(&y[i], vy);
            _mm256_storeu_pd(&r[i], vr);

            #ifdef __FMA__
            vdot0 = _mm256_fmadd_pd(vr, vr, vdot0);
            #else
            vdot0 = _mm256_add_pd(vdot0, _mm256_mul_pd(vr, vr));
            #endif
        }

        __m256d vtotal = _mm256_add_pd(vdot0, vdot1);
        dot = ctools_hsum_avx2(vtotal);
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t valpha = vdupq_n_f64(alpha);
        float64x2_t vdot0 = vdupq_n_f64(0.0);
        float64x2_t vdot1 = vdupq_n_f64(0.0);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vu0 = vld1q_f64(&u[i]);
            float64x2_t vr0 = vld1q_f64(&r[i]);
            float64x2_t vv0 = vld1q_f64(&v[i]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);
            float64x2_t vu1 = vld1q_f64(&u[i + 2]);
            float64x2_t vr1 = vld1q_f64(&r[i + 2]);
            float64x2_t vv1 = vld1q_f64(&v[i + 2]);

            /* y -= alpha * u (fnma: result = acc - a * b) */
            vy0 = vfmsq_f64(vy0, valpha, vu0);
            vy1 = vfmsq_f64(vy1, valpha, vu1);

            /* r -= alpha * v */
            vr0 = vfmsq_f64(vr0, valpha, vv0);
            vr1 = vfmsq_f64(vr1, valpha, vv1);

            vst1q_f64(&y[i], vy0);
            vst1q_f64(&y[i + 2], vy1);
            vst1q_f64(&r[i], vr0);
            vst1q_f64(&r[i + 2], vr1);

            /* Accumulate r·r */
            vdot0 = vfmaq_f64(vdot0, vr0, vr0);
            vdot1 = vfmaq_f64(vdot1, vr1, vr1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vy = vld1q_f64(&y[i]);
            float64x2_t vu = vld1q_f64(&u[i]);
            float64x2_t vr = vld1q_f64(&r[i]);
            float64x2_t vv = vld1q_f64(&v[i]);

            vy = vfmsq_f64(vy, valpha, vu);
            vr = vfmsq_f64(vr, valpha, vv);

            vst1q_f64(&y[i], vy);
            vst1q_f64(&r[i], vr);

            vdot0 = vfmaq_f64(vdot0, vr, vr);
        }

        float64x2_t vtotal = vaddq_f64(vdot0, vdot1);
        dot = ctools_hsum_neon(vtotal);
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        y[i] -= alpha * u[i];
        r[i] -= alpha * v[i];
        dot += r[i] * r[i];
    }

    return dot;
}

/*
 * ctools_simd_fused_axpy2_wdot - Fused weighted variant
 *
 * Performs y -= alpha*u, r -= alpha*v, returns sum(w[i] * r[i] * r[i]).
 */
static inline double ctools_simd_fused_axpy2_wdot(double *y, const double *u,
                                                    double *r, const double *v,
                                                    const double *w,
                                                    double alpha, size_t n)
{
    double dot = 0.0;
    size_t i = 0;

#if CTOOLS_HAS_AVX2
    if (n >= CTOOLS_SIMD_THRESHOLD) {
        __m256d valpha = _mm256_set1_pd(alpha);
        __m256d vdot0 = _mm256_setzero_pd();
        __m256d vdot1 = _mm256_setzero_pd();

        for (; i + 8 <= n; i += 8) {
            __m256d vy0 = _mm256_loadu_pd(&y[i]);
            __m256d vu0 = _mm256_loadu_pd(&u[i]);
            __m256d vr0 = _mm256_loadu_pd(&r[i]);
            __m256d vv0 = _mm256_loadu_pd(&v[i]);
            __m256d vw0 = _mm256_loadu_pd(&w[i]);
            __m256d vy1 = _mm256_loadu_pd(&y[i + 4]);
            __m256d vu1 = _mm256_loadu_pd(&u[i + 4]);
            __m256d vr1 = _mm256_loadu_pd(&r[i + 4]);
            __m256d vv1 = _mm256_loadu_pd(&v[i + 4]);
            __m256d vw1 = _mm256_loadu_pd(&w[i + 4]);

            #ifdef __FMA__
            vy0 = _mm256_fnmadd_pd(valpha, vu0, vy0);
            vy1 = _mm256_fnmadd_pd(valpha, vu1, vy1);
            vr0 = _mm256_fnmadd_pd(valpha, vv0, vr0);
            vr1 = _mm256_fnmadd_pd(valpha, vv1, vr1);
            #else
            vy0 = _mm256_sub_pd(vy0, _mm256_mul_pd(valpha, vu0));
            vy1 = _mm256_sub_pd(vy1, _mm256_mul_pd(valpha, vu1));
            vr0 = _mm256_sub_pd(vr0, _mm256_mul_pd(valpha, vv0));
            vr1 = _mm256_sub_pd(vr1, _mm256_mul_pd(valpha, vv1));
            #endif

            _mm256_storeu_pd(&y[i], vy0);
            _mm256_storeu_pd(&y[i + 4], vy1);
            _mm256_storeu_pd(&r[i], vr0);
            _mm256_storeu_pd(&r[i + 4], vr1);

            /* Accumulate w * r * r */
            __m256d vwr0 = _mm256_mul_pd(vw0, vr0);
            __m256d vwr1 = _mm256_mul_pd(vw1, vr1);
            #ifdef __FMA__
            vdot0 = _mm256_fmadd_pd(vwr0, vr0, vdot0);
            vdot1 = _mm256_fmadd_pd(vwr1, vr1, vdot1);
            #else
            vdot0 = _mm256_add_pd(vdot0, _mm256_mul_pd(vwr0, vr0));
            vdot1 = _mm256_add_pd(vdot1, _mm256_mul_pd(vwr1, vr1));
            #endif
        }

        for (; i + 4 <= n; i += 4) {
            __m256d vy = _mm256_loadu_pd(&y[i]);
            __m256d vu = _mm256_loadu_pd(&u[i]);
            __m256d vr = _mm256_loadu_pd(&r[i]);
            __m256d vv = _mm256_loadu_pd(&v[i]);
            __m256d vw = _mm256_loadu_pd(&w[i]);

            #ifdef __FMA__
            vy = _mm256_fnmadd_pd(valpha, vu, vy);
            vr = _mm256_fnmadd_pd(valpha, vv, vr);
            #else
            vy = _mm256_sub_pd(vy, _mm256_mul_pd(valpha, vu));
            vr = _mm256_sub_pd(vr, _mm256_mul_pd(valpha, vv));
            #endif

            _mm256_storeu_pd(&y[i], vy);
            _mm256_storeu_pd(&r[i], vr);

            __m256d vwr = _mm256_mul_pd(vw, vr);
            #ifdef __FMA__
            vdot0 = _mm256_fmadd_pd(vwr, vr, vdot0);
            #else
            vdot0 = _mm256_add_pd(vdot0, _mm256_mul_pd(vwr, vr));
            #endif
        }

        __m256d vtotal = _mm256_add_pd(vdot0, vdot1);
        dot = ctools_hsum_avx2(vtotal);
    }
#elif CTOOLS_HAS_NEON
    if (n >= 4) {
        float64x2_t valpha = vdupq_n_f64(alpha);
        float64x2_t vdot0 = vdupq_n_f64(0.0);
        float64x2_t vdot1 = vdupq_n_f64(0.0);

        for (; i + 4 <= n; i += 4) {
            float64x2_t vy0 = vld1q_f64(&y[i]);
            float64x2_t vu0 = vld1q_f64(&u[i]);
            float64x2_t vr0 = vld1q_f64(&r[i]);
            float64x2_t vv0 = vld1q_f64(&v[i]);
            float64x2_t vw0 = vld1q_f64(&w[i]);
            float64x2_t vy1 = vld1q_f64(&y[i + 2]);
            float64x2_t vu1 = vld1q_f64(&u[i + 2]);
            float64x2_t vr1 = vld1q_f64(&r[i + 2]);
            float64x2_t vv1 = vld1q_f64(&v[i + 2]);
            float64x2_t vw1 = vld1q_f64(&w[i + 2]);

            vy0 = vfmsq_f64(vy0, valpha, vu0);
            vy1 = vfmsq_f64(vy1, valpha, vu1);
            vr0 = vfmsq_f64(vr0, valpha, vv0);
            vr1 = vfmsq_f64(vr1, valpha, vv1);

            vst1q_f64(&y[i], vy0);
            vst1q_f64(&y[i + 2], vy1);
            vst1q_f64(&r[i], vr0);
            vst1q_f64(&r[i + 2], vr1);

            float64x2_t vwr0 = vmulq_f64(vw0, vr0);
            float64x2_t vwr1 = vmulq_f64(vw1, vr1);
            vdot0 = vfmaq_f64(vdot0, vwr0, vr0);
            vdot1 = vfmaq_f64(vdot1, vwr1, vr1);
        }

        for (; i + 2 <= n; i += 2) {
            float64x2_t vy = vld1q_f64(&y[i]);
            float64x2_t vu = vld1q_f64(&u[i]);
            float64x2_t vr = vld1q_f64(&r[i]);
            float64x2_t vv = vld1q_f64(&v[i]);
            float64x2_t vw = vld1q_f64(&w[i]);

            vy = vfmsq_f64(vy, valpha, vu);
            vr = vfmsq_f64(vr, valpha, vv);

            vst1q_f64(&y[i], vy);
            vst1q_f64(&r[i], vr);

            float64x2_t vwr = vmulq_f64(vw, vr);
            vdot0 = vfmaq_f64(vdot0, vwr, vr);
        }

        float64x2_t vtotal = vaddq_f64(vdot0, vdot1);
        dot = ctools_hsum_neon(vtotal);
    }
#endif

    /* Scalar tail */
    for (; i < n; i++) {
        y[i] -= alpha * u[i];
        r[i] -= alpha * v[i];
        dot += w[i] * r[i] * r[i];
    }

    return dot;
}

#endif /* CTOOLS_SIMD_H */

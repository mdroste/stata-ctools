/*
 * ctools_unroll.h
 *
 * Loop unrolling abstraction for ctools.
 * Provides macros and inline functions for K-way unrolled reductions.
 *
 * Usage:
 *   - Change CTOOLS_UNROLL_K to experiment with different unroll factors
 *   - Use CTOOLS_DOT_PRODUCT() for standard dot products
 *   - Use CTOOLS_REDUCE_SUM() for custom reduction operations
 *
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_UNROLL_H
#define CTOOLS_UNROLL_H

#include "stplugin.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

/*
 * CTOOLS_UNROLL_K: Number of independent accumulator chains.
 *
 * Supported values: 4, 8 (default), 16
 *
 * Trade-offs:
 *   K=4:  Lower register pressure, works well on older CPUs
 *   K=8:  Good balance for modern out-of-order CPUs (default)
 *   K=16: Maximum ILP, but may spill registers on some architectures
 *
 * To experiment, redefine before including this header:
 *   #define CTOOLS_UNROLL_K 4
 *   #include "ctools_unroll.h"
 */
#ifndef CTOOLS_UNROLL_K
#define CTOOLS_UNROLL_K 8
#endif

/* Compile-time validation */
#if CTOOLS_UNROLL_K != 4 && CTOOLS_UNROLL_K != 8 && CTOOLS_UNROLL_K != 16
#error "CTOOLS_UNROLL_K must be 4, 8, or 16"
#endif

/* ============================================================================
 * Internal Macros - Do Not Use Directly
 * ============================================================================ */

/* Declare K accumulators initialized to zero */
#if CTOOLS_UNROLL_K == 4
#define CTOOLS_DECLARE_ACCUM(prefix) \
    ST_double prefix##0 = 0.0, prefix##1 = 0.0, \
              prefix##2 = 0.0, prefix##3 = 0.0
#elif CTOOLS_UNROLL_K == 8
#define CTOOLS_DECLARE_ACCUM(prefix) \
    ST_double prefix##0 = 0.0, prefix##1 = 0.0, \
              prefix##2 = 0.0, prefix##3 = 0.0, \
              prefix##4 = 0.0, prefix##5 = 0.0, \
              prefix##6 = 0.0, prefix##7 = 0.0
#elif CTOOLS_UNROLL_K == 16
#define CTOOLS_DECLARE_ACCUM(prefix) \
    ST_double prefix##0 = 0.0, prefix##1 = 0.0, \
              prefix##2 = 0.0, prefix##3 = 0.0, \
              prefix##4 = 0.0, prefix##5 = 0.0, \
              prefix##6 = 0.0, prefix##7 = 0.0, \
              prefix##8 = 0.0, prefix##9 = 0.0, \
              prefix##10 = 0.0, prefix##11 = 0.0, \
              prefix##12 = 0.0, prefix##13 = 0.0, \
              prefix##14 = 0.0, prefix##15 = 0.0
#endif

/* K-way unrolled accumulation: prefix_k += x[i+k] * y[i+k] */
#if CTOOLS_UNROLL_K == 4
#define CTOOLS_ACCUM_DOT(prefix, x, y, i) \
    do { \
        prefix##0 += (x)[(i)]     * (y)[(i)];     \
        prefix##1 += (x)[(i) + 1] * (y)[(i) + 1]; \
        prefix##2 += (x)[(i) + 2] * (y)[(i) + 2]; \
        prefix##3 += (x)[(i) + 3] * (y)[(i) + 3]; \
    } while (0)
#elif CTOOLS_UNROLL_K == 8
#define CTOOLS_ACCUM_DOT(prefix, x, y, i) \
    do { \
        prefix##0 += (x)[(i)]     * (y)[(i)];     \
        prefix##1 += (x)[(i) + 1] * (y)[(i) + 1]; \
        prefix##2 += (x)[(i) + 2] * (y)[(i) + 2]; \
        prefix##3 += (x)[(i) + 3] * (y)[(i) + 3]; \
        prefix##4 += (x)[(i) + 4] * (y)[(i) + 4]; \
        prefix##5 += (x)[(i) + 5] * (y)[(i) + 5]; \
        prefix##6 += (x)[(i) + 6] * (y)[(i) + 6]; \
        prefix##7 += (x)[(i) + 7] * (y)[(i) + 7]; \
    } while (0)
#elif CTOOLS_UNROLL_K == 16
#define CTOOLS_ACCUM_DOT(prefix, x, y, i) \
    do { \
        prefix##0  += (x)[(i)]      * (y)[(i)];      \
        prefix##1  += (x)[(i) + 1]  * (y)[(i) + 1];  \
        prefix##2  += (x)[(i) + 2]  * (y)[(i) + 2];  \
        prefix##3  += (x)[(i) + 3]  * (y)[(i) + 3];  \
        prefix##4  += (x)[(i) + 4]  * (y)[(i) + 4];  \
        prefix##5  += (x)[(i) + 5]  * (y)[(i) + 5];  \
        prefix##6  += (x)[(i) + 6]  * (y)[(i) + 6];  \
        prefix##7  += (x)[(i) + 7]  * (y)[(i) + 7];  \
        prefix##8  += (x)[(i) + 8]  * (y)[(i) + 8];  \
        prefix##9  += (x)[(i) + 9]  * (y)[(i) + 9];  \
        prefix##10 += (x)[(i) + 10] * (y)[(i) + 10]; \
        prefix##11 += (x)[(i) + 11] * (y)[(i) + 11]; \
        prefix##12 += (x)[(i) + 12] * (y)[(i) + 12]; \
        prefix##13 += (x)[(i) + 13] * (y)[(i) + 13]; \
        prefix##14 += (x)[(i) + 14] * (y)[(i) + 14]; \
        prefix##15 += (x)[(i) + 15] * (y)[(i) + 15]; \
    } while (0)
#endif

/* K-way unrolled self dot: prefix_k += x[i+k] * x[i+k] */
#if CTOOLS_UNROLL_K == 4
#define CTOOLS_ACCUM_DOT_SELF(prefix, x, i) \
    do { \
        prefix##0 += (x)[(i)]     * (x)[(i)];     \
        prefix##1 += (x)[(i) + 1] * (x)[(i) + 1]; \
        prefix##2 += (x)[(i) + 2] * (x)[(i) + 2]; \
        prefix##3 += (x)[(i) + 3] * (x)[(i) + 3]; \
    } while (0)
#elif CTOOLS_UNROLL_K == 8
#define CTOOLS_ACCUM_DOT_SELF(prefix, x, i) \
    do { \
        prefix##0 += (x)[(i)]     * (x)[(i)];     \
        prefix##1 += (x)[(i) + 1] * (x)[(i) + 1]; \
        prefix##2 += (x)[(i) + 2] * (x)[(i) + 2]; \
        prefix##3 += (x)[(i) + 3] * (x)[(i) + 3]; \
        prefix##4 += (x)[(i) + 4] * (x)[(i) + 4]; \
        prefix##5 += (x)[(i) + 5] * (x)[(i) + 5]; \
        prefix##6 += (x)[(i) + 6] * (x)[(i) + 6]; \
        prefix##7 += (x)[(i) + 7] * (x)[(i) + 7]; \
    } while (0)
#elif CTOOLS_UNROLL_K == 16
#define CTOOLS_ACCUM_DOT_SELF(prefix, x, i) \
    do { \
        prefix##0  += (x)[(i)]      * (x)[(i)];      \
        prefix##1  += (x)[(i) + 1]  * (x)[(i) + 1];  \
        prefix##2  += (x)[(i) + 2]  * (x)[(i) + 2];  \
        prefix##3  += (x)[(i) + 3]  * (x)[(i) + 3];  \
        prefix##4  += (x)[(i) + 4]  * (x)[(i) + 4];  \
        prefix##5  += (x)[(i) + 5]  * (x)[(i) + 5];  \
        prefix##6  += (x)[(i) + 6]  * (x)[(i) + 6];  \
        prefix##7  += (x)[(i) + 7]  * (x)[(i) + 7];  \
        prefix##8  += (x)[(i) + 8]  * (x)[(i) + 8];  \
        prefix##9  += (x)[(i) + 9]  * (x)[(i) + 9];  \
        prefix##10 += (x)[(i) + 10] * (x)[(i) + 10]; \
        prefix##11 += (x)[(i) + 11] * (x)[(i) + 11]; \
        prefix##12 += (x)[(i) + 12] * (x)[(i) + 12]; \
        prefix##13 += (x)[(i) + 13] * (x)[(i) + 13]; \
        prefix##14 += (x)[(i) + 14] * (x)[(i) + 14]; \
        prefix##15 += (x)[(i) + 15] * (x)[(i) + 15]; \
    } while (0)
#endif

/*
 * Tree reduction for numerical stability.
 * Pairs accumulators to minimize floating-point error propagation.
 */
#if CTOOLS_UNROLL_K == 4
#define CTOOLS_TREE_REDUCE(prefix) \
    (((prefix##0) + (prefix##2)) + ((prefix##1) + (prefix##3)))
#elif CTOOLS_UNROLL_K == 8
#define CTOOLS_TREE_REDUCE(prefix) \
    ((((prefix##0) + (prefix##4)) + ((prefix##1) + (prefix##5))) + \
     (((prefix##2) + (prefix##6)) + ((prefix##3) + (prefix##7))))
#elif CTOOLS_UNROLL_K == 16
#define CTOOLS_TREE_REDUCE(prefix) \
    (((((prefix##0) + (prefix##8))  + ((prefix##1) + (prefix##9)))  + \
      (((prefix##2) + (prefix##10)) + ((prefix##3) + (prefix##11)))) + \
     ((((prefix##4) + (prefix##12)) + ((prefix##5) + (prefix##13))) + \
      (((prefix##6) + (prefix##14)) + ((prefix##7) + (prefix##15)))))
#endif

/* ============================================================================
 * Public API - Inline Functions
 * ============================================================================ */

/*
 * ctools_dot_unrolled - K-way unrolled dot product
 *
 * Computes x . y using CTOOLS_UNROLL_K independent accumulator chains
 * with tree reduction for numerical stability.
 *
 * Parameters:
 *   x, y - Input arrays (must have at least N elements)
 *   N    - Number of elements
 *
 * Returns:
 *   The dot product sum(x[i] * y[i]) for i in [0, N)
 */
static inline ST_double ctools_dot_unrolled(
    const ST_double * restrict x,
    const ST_double * restrict y,
    ST_int N)
{
    ST_int i;
    CTOOLS_DECLARE_ACCUM(sum);
    ST_int N_main = N - (N % CTOOLS_UNROLL_K);

    /* K-way unrolled main loop */
    for (i = 0; i < N_main; i += CTOOLS_UNROLL_K) {
        CTOOLS_ACCUM_DOT(sum, x, y, i);
    }

    /* Scalar remainder */
    for (; i < N; i++) {
        sum0 += x[i] * y[i];
    }

    return CTOOLS_TREE_REDUCE(sum);
}

/*
 * ctools_dot_self_unrolled - K-way unrolled self dot product (norm squared)
 *
 * Computes ||x||^2 = x . x using CTOOLS_UNROLL_K independent accumulator chains.
 *
 * Parameters:
 *   x - Input array (must have at least N elements)
 *   N - Number of elements
 *
 * Returns:
 *   The squared norm sum(x[i] * x[i]) for i in [0, N)
 */
static inline ST_double ctools_dot_self_unrolled(
    const ST_double * restrict x,
    ST_int N)
{
    ST_int i;
    CTOOLS_DECLARE_ACCUM(sum);
    ST_int N_main = N - (N % CTOOLS_UNROLL_K);

    /* K-way unrolled main loop */
    for (i = 0; i < N_main; i += CTOOLS_UNROLL_K) {
        CTOOLS_ACCUM_DOT_SELF(sum, x, i);
    }

    /* Scalar remainder */
    for (; i < N; i++) {
        sum0 += x[i] * x[i];
    }

    return CTOOLS_TREE_REDUCE(sum);
}

/* ============================================================================
 * Convenience Macros for Custom Reductions
 * ============================================================================
 *
 * For operations beyond dot products, use these building blocks:
 *
 * Example - weighted dot product:
 *
 *   ST_double weighted_dot(const ST_double *x, const ST_double *y,
 *                          const ST_double *w, ST_int N) {
 *       ST_int i;
 *       CTOOLS_DECLARE_ACCUM(sum);
 *       ST_int N_main = N - (N % CTOOLS_UNROLL_K);
 *
 *       for (i = 0; i < N_main; i += CTOOLS_UNROLL_K) {
 *           // Custom K-way accumulation (must match CTOOLS_UNROLL_K)
 *           #if CTOOLS_UNROLL_K >= 4
 *           sum0 += w[i]     * x[i]     * y[i];
 *           sum1 += w[i + 1] * x[i + 1] * y[i + 1];
 *           sum2 += w[i + 2] * x[i + 2] * y[i + 2];
 *           sum3 += w[i + 3] * x[i + 3] * y[i + 3];
 *           #endif
 *           #if CTOOLS_UNROLL_K >= 8
 *           sum4 += w[i + 4] * x[i + 4] * y[i + 4];
 *           sum5 += w[i + 5] * x[i + 5] * y[i + 5];
 *           sum6 += w[i + 6] * x[i + 6] * y[i + 6];
 *           sum7 += w[i + 7] * x[i + 7] * y[i + 7];
 *           #endif
 *           // ... etc for K=16
 *       }
 *       for (; i < N; i++) sum0 += w[i] * x[i] * y[i];
 *       return CTOOLS_TREE_REDUCE(sum);
 *   }
 */

#endif /* CTOOLS_UNROLL_H */

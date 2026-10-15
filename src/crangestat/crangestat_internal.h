#ifndef CRANGESTAT_INTERNAL_H
#define CRANGESTAT_INTERNAL_H
/*
    crangestat_impl.c
    High-performance range statistics for Stata

    Performance optimizations:
    - Parallel IPS4o sort (via ctools infrastructure)
    - Binary search for O(log n) window finding
    - Run-length optimization: same-key observations share results
    - Single-pass statistics without data copying
    - SIMD-accelerated window statistics
    - Small group batching to reduce parallel overhead
    - OpenMP parallelization across observations/runs

    Author: ctools project
    License: MIT
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Checked SPI writes resolve positions in the plugin varlist, even after
 * estimation has introduced internal Stata variables. */
#define SD_SAFEMODE
#include "stplugin.h"
#include "ctools_simd.h"
#include "ctools_types.h"
#include "ctools_parse.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_order.h"
#include "ctools_select.h"

#define CRANGESTAT_MODULE "crangestat"

/* Thresholds for optimizations */
#define CRANGESTAT_RUN_OPT_THRESHOLD 1000    /* Min group size for run optimization */
#define CRANGESTAT_BATCH_THRESHOLD 100       /* Max group size for batching */
#define CRANGESTAT_PARALLEL_THRESHOLD 10000  /* Min work for parallelization */
#define CRANGESTAT_PARALLEL_GROUPS_MIN 8     /* Min groups for cross-group parallelism */
#define CRANGESTAT_SMALL_GROUP_THRESHOLD 500 /* Max avg group size for cross-group parallelism */
#define CRANGESTAT_PREFIX_SUM_THRESHOLD 64   /* Min window size to use prefix sums */
#define CRANGESTAT_SPARSE_TABLE_LOG 20       /* Max log2(n) for sparse table */

/* Statistics identifiers */
typedef enum
{
    STAT_COUNT = 0,
    STAT_MEAN,
    STAT_SUM,
    STAT_MIN,
    STAT_MAX,
    STAT_SD,
    STAT_VARIANCE,
    STAT_MEDIAN,
    STAT_IQR,
    STAT_FIRST,
    STAT_LAST,
    STAT_FIRSTNM,
    STAT_LASTNM,
    STAT_P1,
    STAT_P5,
    STAT_P10,
    STAT_P25,
    STAT_P75,
    STAT_P90,
    STAT_P95,
    STAT_P99,
    STAT_SKEWNESS,
    STAT_KURTOSIS,
    STAT_COUNT_ALL /* Total number of statistic types */
} stat_type;

/* Configuration */
typedef struct
{
    int exclude_self;
    int verbose;
    double interval_low;  /* Offset for low bound (can be . for -inf) */
    double interval_high; /* Offset for high bound (can be . for +inf) */
    int low_is_missing;   /* 1 if low bound is . (meaning -infinity) */
    int high_is_missing;  /* 1 if high bound is . (meaning +infinity) */
} crangestat_config;

/* Statistic specification */
typedef struct
{
    stat_type type;
    int source_var_idx; /* Index into loaded source variables */
    int result_var_idx; /* 1-based Stata variable index for result */
} stat_spec;

/* ===========================================================================
   Binary search for window boundaries
   =========================================================================== */

/* Find first index where key[i] >= target */

/* Find first index where key[i] > target */

/* ===========================================================================
   Quickselect for median/percentiles (shared ctools_quickselect_double,
   which stays O(n) on tied windows)
   =========================================================================== */

/* ===========================================================================
   Group detection
   =========================================================================== */

typedef struct
{
    size_t start;
    size_t count;
} group_info;

/* ===========================================================================
   Prefix Sum Arrays for O(1) Range Queries

   For sum, count, mean, sd, variance, we precompute prefix sums:
   - prefix_sum[i] = sum of data[0..i-1] (exclusive end)
   - prefix_sum2[i] = sum of data[0..i-1]^2
   - prefix_count[i] = count of non-missing in [0..i-1]

   Range query: sum[a..b) = prefix[b] - prefix[a]
   With excludeself: subtract self contribution
   =========================================================================== */

typedef struct
{
    double *prefix_sum;   /* Prefix sums: prefix_sum[i] = sum of [0,i) */
    double *prefix_sum2;  /* Prefix sum of squares */
    size_t *prefix_count; /* Prefix count of non-missing */
    size_t nobs;
    double origin; /* Center moments before accumulating squares. */
} prefix_arrays;

/* Build prefix arrays for a source variable within a group */

/* Query sum for range [win_start, win_end) within group starting at g_start */

/* ===========================================================================
   Sparse Table for O(1) Range Min/Max Queries

   sparse_min[i][j] = min of range [i, i + 2^j)
   sparse_max[i][j] = max of range [i, i + 2^j)

   Query [l, r]: k = log2(r - l), answer = min(sparse[l][k], sparse[r - 2^k][k])
   =========================================================================== */

typedef struct
{
    double **sparse_min; /* sparse_min[i][j] = min of [i, i + 2^j) */
    double **sparse_max; /* sparse_max[i][j] = max of [i, i + 2^j) */
    int *log_table;      /* Precomputed log2 values */
    size_t nobs;
    int max_log;
} sparse_table;

/* Build sparse table for range min/max queries */

/* Query range min for [l, r) - l and r are group-relative.
 * Levels stop at max_log to bound memory, so a range longer than
 * 2^(max_log+1) is covered by several top-level blocks; shorter ranges use
 * the usual two overlapping blocks. */

/* Query range max for [l, r) - l and r are group-relative (see query_range_min) */

/* Free sparse table */

/* Allocate sparse table buffers without filling data.
 * max_count is the maximum group size; actual groups may be smaller.
 * Returns 0 on success, -1 on allocation failure. */

/* Rebuild sparse table data for a specific group using pre-allocated buffers.
 * st must have been allocated with alloc_sparse_table(st, max_count) where
 * max_count >= count. Pointer stride remains max_count (the allocation capacity).
 * Only count entries per level are meaningful; queries must use indices < count. */

/* ===========================================================================
   Statistics computation for a window
   =========================================================================== */

/* Check if stat type needs work array (percentiles, etc.) */

/* Check if stat can use O(1) prefix sum queries */

/* Check if stat can use O(1) sparse table queries */

/* Two-pass centered moments for small windows and ill-conditioned prefix
 * differences. Long double helps on platforms where it exceeds double precision. */

/*
    compute_stat_with_prefix - O(1) stat computation using prefix sums

    For count, sum, mean, sd, variance, this is O(1) instead of O(window_size)!
    This is a massive speedup for large windows.
*/

/*
    compute_minmax_with_sparse - O(1) min/max using sparse table

    For min/max without excludeself, this is O(1) instead of O(window_size)!
    With excludeself, we check if result equals self and rescan if needed.
*/

/* SIMD-accelerated sum for window statistics */

/* Compute simple statistics directly without copying to work array */

/* Compute statistics that need work array (percentiles, etc.) */

/* ===========================================================================
   Argument parsing
   =========================================================================== */

/* ===========================================================================
   Main Entry Point
   =========================================================================== */

/* Private kernels: caller owns samples, workspaces and output buffers.
 * Matching/window traversal preserve stable order and existing definitions. */
size_t crangestat_lower_bound(const double *key, size_t start, size_t end, double target);
size_t crangestat_upper_bound(const double *key, size_t start, size_t end, double target);
double crangestat_compute_percentile(double *work, size_t n, double pctl);
void crangestat_build_prefix_arrays(const double *data, size_t start, size_t count,
                                    prefix_arrays *pa);
void crangestat_query_prefix_sums(const prefix_arrays *pa, size_t g_start, size_t win_start,
                                  size_t win_end, double *out_sum, double *out_sum2,
                                  size_t *out_count);
int crangestat_build_sparse_table(const double *data, size_t start, size_t count, sparse_table *st,
                                  const double miss);
double crangestat_query_range_min(const sparse_table *st, size_t l, size_t r);
double crangestat_query_range_max(const sparse_table *st, size_t l, size_t r);
void crangestat_free_sparse_table(sparse_table *st);
int crangestat_alloc_sparse_table(sparse_table *st, size_t max_count);
void crangestat_rebuild_sparse_table(sparse_table *st, const double *data, size_t start,
                                     size_t count, double miss);
int crangestat_stat_needs_work_array(stat_type type);
int crangestat_stat_uses_prefix_sums(stat_type type);
int crangestat_stat_uses_sparse_table(stat_type type);
double crangestat_stable_window_variance(const double *data, size_t start, size_t end, size_t self,
                                         int exclude_self);
double crangestat_compute_stat_with_prefix(stat_type type, const prefix_arrays *pa,
                                           const double *data, size_t g_start, size_t win_start,
                                           size_t win_end, size_t self_idx, int exclude_self);
double crangestat_compute_minmax_with_sparse(stat_type type, const sparse_table *st,
                                             const double *data, size_t g_start, size_t win_start,
                                             size_t win_end, size_t self_idx, int exclude_self);
void crangestat_simd_sum_minmax(const double *data, size_t start, size_t end, double *out_sum,
                                double *out_sum2, double *out_min, double *out_max, size_t *out_n);
double crangestat_compute_window_stat_simple(stat_type type, const double *data, size_t win_start,
                                             size_t win_end, size_t self_idx, int exclude_self);
double crangestat_compute_window_stat_complex(stat_type type, const double *data, size_t win_start,
                                              size_t win_end, size_t self_idx, int exclude_self,
                                              double *work);
#endif

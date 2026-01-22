/*
 * ctools_select.h
 * Unified quickselect algorithms for O(n) k-th element selection
 *
 * Provides:
 * - ctools_quickselect_double: Select k-th smallest from double array
 * - ctools_quickselect_indexed: Select k-th smallest via index array
 *
 * All functions use median-of-three pivot selection for better performance
 * on partially sorted data.
 */

#ifndef CTOOLS_SELECT_H
#define CTOOLS_SELECT_H

#include <stddef.h>
#include "stplugin.h"

/* ============================================================================
 * Direct quickselect for double arrays
 *
 * Finds the k-th smallest element in O(n) expected time.
 * Modifies array in place (partial sorting).
 * Uses insertion sort for small subarrays (< 10 elements).
 *
 * @param arr   Array of doubles (modified in place)
 * @param n     Number of elements
 * @param k     Index of element to find (0-based)
 * @return      Value of k-th smallest element, or SV_missval if n == 0
 * ============================================================================ */

static inline void ctools_swap_double(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline size_t ctools_partition_double(double *arr, size_t left, size_t right)
{
    size_t mid = left + (right - left) / 2;

    /* Median-of-three pivot selection */
    if (arr[mid] < arr[left]) ctools_swap_double(&arr[left], &arr[mid]);
    if (arr[right] < arr[left]) ctools_swap_double(&arr[left], &arr[right]);
    if (arr[right] < arr[mid]) ctools_swap_double(&arr[mid], &arr[right]);

    /* Move pivot to right-1 */
    ctools_swap_double(&arr[mid], &arr[right - 1]);
    double pivot = arr[right - 1];

    size_t i = left;
    for (size_t j = left; j < right - 1; j++) {
        if (arr[j] <= pivot) {
            ctools_swap_double(&arr[i], &arr[j]);
            i++;
        }
    }
    ctools_swap_double(&arr[i], &arr[right - 1]);
    return i;
}

static inline double ctools_quickselect_double(double *arr, size_t n, size_t k)
{
    if (n == 0) return SV_missval;
    if (n == 1) return arr[0];
    if (k >= n) k = n - 1;

    size_t left = 0;
    size_t right = n - 1;

    while (left < right) {
        /* Insertion sort for small subarrays */
        if (right - left < 10) {
            for (size_t i = left + 1; i <= right; i++) {
                double key = arr[i];
                size_t j = i;
                while (j > left && arr[j - 1] > key) {
                    arr[j] = arr[j - 1];
                    j--;
                }
                arr[j] = key;
            }
            return arr[k];
        }

        size_t pivot_idx = ctools_partition_double(arr, left, right);
        if (k == pivot_idx) return arr[k];
        else if (k < pivot_idx) right = pivot_idx - 1;
        else left = pivot_idx + 1;
    }
    return arr[left];
}

/* ============================================================================
 * Quickselect with Stata types (ST_double, ST_int)
 *
 * Same algorithm but using Stata's type definitions for consistency
 * with the rest of ctools.
 * ============================================================================ */

static inline void ctools_quickselect_st(ST_double *arr, ST_int lo, ST_int hi, ST_int k)
{
    while (lo < hi) {
        ST_int mid = lo + (hi - lo) / 2;

        /* Median-of-three pivot selection */
        if (arr[mid] < arr[lo]) {
            ST_double tmp = arr[lo]; arr[lo] = arr[mid]; arr[mid] = tmp;
        }
        if (arr[hi] < arr[lo]) {
            ST_double tmp = arr[lo]; arr[lo] = arr[hi]; arr[hi] = tmp;
        }
        if (arr[mid] < arr[hi]) {
            ST_double tmp = arr[mid]; arr[mid] = arr[hi]; arr[hi] = tmp;
        }

        /* Pivot is now at hi */
        ST_double pivot = arr[hi];
        ST_int i = lo - 1;

        for (ST_int j = lo; j < hi; j++) {
            if (arr[j] <= pivot) {
                i++;
                ST_double tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
            }
        }

        ST_double tmp = arr[i + 1]; arr[i + 1] = arr[hi]; arr[hi] = tmp;
        ST_int pivot_idx = i + 1;

        if (pivot_idx == k) {
            return;
        } else if (k < pivot_idx) {
            hi = pivot_idx - 1;
        } else {
            lo = pivot_idx + 1;
        }
    }
}

/* ============================================================================
 * Indexed quickselect
 *
 * Reorders an index array such that idx[k] points to the k-th smallest
 * value in vals. Does not modify vals array.
 *
 * @param vals  Array of values (read-only)
 * @param idx   Array of indices into vals (modified in place)
 * @param lo    Start index in idx array
 * @param hi    End index in idx array (inclusive)
 * @return      Final position of pivot (the k-th element position)
 * ============================================================================ */

static inline void ctools_swap_int(ST_int *a, ST_int *b)
{
    ST_int t = *a; *a = *b; *b = t;
}

static inline ST_int ctools_partition_indexed(const ST_double *vals, ST_int *idx,
                                               ST_int lo, ST_int hi)
{
    ST_int mid = lo + (hi - lo) / 2;

    /* Median-of-three pivot selection based on values */
    if (vals[idx[mid]] < vals[idx[lo]]) ctools_swap_int(&idx[lo], &idx[mid]);
    if (vals[idx[hi]] < vals[idx[lo]]) ctools_swap_int(&idx[lo], &idx[hi]);
    if (vals[idx[hi]] < vals[idx[mid]]) ctools_swap_int(&idx[mid], &idx[hi]);

    /* Move pivot to hi-1 */
    ctools_swap_int(&idx[mid], &idx[hi - 1]);
    ST_double pivot = vals[idx[hi - 1]];

    ST_int i = lo;
    ST_int j = hi - 1;

    for (;;) {
        while (vals[idx[++i]] < pivot) {}
        while (vals[idx[--j]] > pivot) { if (j == lo) break; }
        if (i >= j) break;
        ctools_swap_int(&idx[i], &idx[j]);
    }

    ctools_swap_int(&idx[i], &idx[hi - 1]);
    return i;
}

static inline void ctools_quickselect_indexed(const ST_double *vals, ST_int *idx,
                                               ST_int lo, ST_int hi, ST_int k)
{
    while (lo < hi) {
        if (hi - lo < 3) {
            /* Handle tiny arrays directly */
            if (hi - lo == 1) {
                if (vals[idx[lo]] > vals[idx[hi]]) ctools_swap_int(&idx[lo], &idx[hi]);
            } else { /* hi - lo == 2 */
                if (vals[idx[lo]] > vals[idx[lo + 1]]) ctools_swap_int(&idx[lo], &idx[lo + 1]);
                if (vals[idx[lo + 1]] > vals[idx[hi]]) ctools_swap_int(&idx[lo + 1], &idx[hi]);
                if (vals[idx[lo]] > vals[idx[lo + 1]]) ctools_swap_int(&idx[lo], &idx[lo + 1]);
            }
            return;
        }

        ST_int pivot_idx = ctools_partition_indexed(vals, idx, lo, hi);
        if (pivot_idx == k) return;
        else if (k < pivot_idx) hi = pivot_idx - 1;
        else lo = pivot_idx + 1;
    }
}

/* ============================================================================
 * Utility: Find two order statistics efficiently
 *
 * Uses quickselect twice, leveraging partial ordering from first call.
 * ============================================================================ */

static inline void ctools_find_two_percentiles(double *arr, size_t n,
                                                size_t k_lo, size_t k_hi,
                                                double *val_lo, double *val_hi)
{
    if (n == 0) {
        *val_lo = SV_missval;
        *val_hi = SV_missval;
        return;
    }

    if (k_lo >= n) k_lo = n - 1;
    if (k_hi >= n) k_hi = n - 1;

    /* Find lower percentile first */
    *val_lo = ctools_quickselect_double(arr, n, k_lo);

    /* Find upper percentile - array is partially sorted, so search from k_lo */
    if (k_hi > k_lo) {
        *val_hi = ctools_quickselect_double(arr + k_lo, n - k_lo, k_hi - k_lo);
    } else {
        *val_hi = *val_lo;
    }
}

#endif /* CTOOLS_SELECT_H */

/*
 * ctools_select.h
 * Quickselect for O(n) expected k-th element selection
 *
 * Provides ctools_quickselect_double, which uses median-of-three pivot
 * selection for better performance on partially sorted data, a partition
 * that splits runs of equal keys evenly (tied data stays O(n)), and a
 * heapselect fallback that bounds the worst case at O(n log n).
 */

#ifndef CTOOLS_SELECT_H
#define CTOOLS_SELECT_H

#include <stddef.h>
#include <stdint.h>
#include "stplugin.h"

/* Partition work allowed before the heapselect fallback, as a multiple of n
 * (at least 1). Expected work is 2-3n, so only adversarial inputs reach it. */
#ifndef CTOOLS_SELECT_WORK_FACTOR
#define CTOOLS_SELECT_WORK_FACTOR 8
#endif

/* ============================================================================
 * Direct quickselect for double arrays
 *
 * Finds the k-th smallest element in O(n) expected time, O(n log n) worst case.
 * Modifies array in place (partial sorting). On return arr[k] is the k-th
 * smallest value, arr[0..k-1] <= arr[k] <= arr[k+1..n-1]; callers rely on
 * this to take further order statistics from either side.
 * Uses insertion sort for small subarrays (< 10 elements).
 * Values must not be NaN (callers pass nonmissing values).
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

/* Partition arr[left..right] (requires right - left >= 2) around a
 * median-of-three pivot; returns the pivot's final index p, with
 * arr[left..p-1] <= arr[p] <= arr[p+1..right]. */
static inline size_t ctools_partition_double(double *arr, size_t left, size_t right)
{
    size_t mid = left + (right - left) / 2;
    size_t q = (right - left) / 4;

    /* Median-of-three pivot selection over the quartile points and the
     * middle (quartile points moved to the ends first); unlike the end
     * points, these are not fooled by organ-pipe or V-shaped data */
    ctools_swap_double(&arr[left], &arr[left + q]);
    ctools_swap_double(&arr[right], &arr[right - q]);
    if (arr[mid] < arr[left]) ctools_swap_double(&arr[left], &arr[mid]);
    if (arr[right] < arr[left]) ctools_swap_double(&arr[left], &arr[right]);
    if (arr[right] < arr[mid]) ctools_swap_double(&arr[mid], &arr[right]);

    /* Move pivot to right-1; arr[left] <= pivot <= arr[right] stop the scans */
    ctools_swap_double(&arr[mid], &arr[right - 1]);
    double pivot = arr[right - 1];

    /* Hoare scans that stop on keys equal to the pivot, so a run of ties is
     * split evenly between the two sides instead of being peeled off one
     * key per pass (which made tied data quadratic). */
    size_t i = left;
    size_t j = right - 1;
    for (;;) {
        while (arr[++i] < pivot) {}
        while (pivot < arr[--j]) {}
        if (i >= j) break;
        ctools_swap_double(&arr[i], &arr[j]);
    }
    ctools_swap_double(&arr[i], &arr[right - 1]);
    return i;
}

/* Sift h[i] down a max-heap of size elements */
static inline void ctools_sift_down_double(double *h, size_t i, size_t size)
{
    double x = h[i];
    for (;;) {
        size_t c = 2 * i + 1;
        if (c >= size) break;
        if (c + 1 < size && h[c + 1] > h[c]) c++;
        if (!(h[c] > x)) break;
        h[i] = h[c];
        i = c;
    }
    h[i] = x;
}

/* Heapselect on arr[left..right] (left <= k <= right): keep the k-left+1
 * smallest values in a max-heap at arr[left..k], then move its root (the
 * k-th smallest) to arr[k]. O(m log m), same postcondition as quickselect. */
static inline void ctools_heapselect_double(double *arr, size_t left,
                                            size_t right, size_t k)
{
    double *h = arr + left;
    size_t size = k - left + 1;

    for (size_t i = size / 2; i-- > 0; )
        ctools_sift_down_double(h, i, size);
    for (size_t i = k + 1; i <= right; i++) {
        if (arr[i] < h[0]) {
            ctools_swap_double(&arr[i], &h[0]);
            ctools_sift_down_double(h, 0, size);
        }
    }
    ctools_swap_double(&h[0], &arr[k]);
}

static inline double ctools_quickselect_double(double *arr, size_t n, size_t k)
{
    if (n == 0) return SV_missval;
    if (n == 1) return arr[0];
    if (k >= n) k = n - 1;

    size_t left = 0;
    size_t right = n - 1;
    size_t budget = n <= SIZE_MAX / CTOOLS_SELECT_WORK_FACTOR
                  ? n * CTOOLS_SELECT_WORK_FACTOR : SIZE_MAX;

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

        /* Bound the total partition work; if pivots keep failing to shrink
         * the range, finish it with heapselect */
        size_t len = right - left + 1;
        if (len > budget) {
            ctools_heapselect_double(arr, left, right, k);
            return arr[k];
        }
        budget -= len;

        size_t pivot_idx = ctools_partition_double(arr, left, right);
        if (k == pivot_idx) return arr[k];
        else if (k < pivot_idx) right = pivot_idx - 1;
        else left = pivot_idx + 1;
    }
    return arr[left];
}

#endif /* CTOOLS_SELECT_H */

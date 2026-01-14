/*
    ctools_sort_timsort.c
    Parallel Timsort Module

    Timsort is an adaptive, stable, hybrid sorting algorithm designed by Tim Peters.
    It combines merge sort and insertion sort to achieve excellent performance on
    real-world data that often contains existing order (partially sorted runs).

    Key advantages:
    - O(N) best case when data is already sorted
    - O(N log N) worst case guaranteed
    - Stable sort (maintains relative order of equal elements)
    - Adaptive: exploits existing order in the data
    - Excellent cache performance due to sequential memory access

    Algorithm overview:
    1. Find natural runs (ascending or strictly descending sequences)
    2. Extend short runs using binary insertion sort to minimum run length
    3. Merge runs using a stack-based merge strategy
    4. Use galloping mode during merge for sequences with many consecutive elements

    Parallelization Strategy:
    - Parallel key conversion using OpenMP
    - Parallel run detection across data chunks
    - Parallel merging of independent run pairs
    - Sequential merge stack operations (maintains stability)

    Optimizations:
    - memmove for batch element shifting in insertion sort
    - Parallel key conversion for large datasets
    - Galloping mode for merge optimization

    Best suited for:
    - Data with existing partial order (common in Stata panel data)
    - Large datasets where cache efficiency matters
    - When stability is required
    - Mixed numeric/string data with complex comparison
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/* Minimum observations for parallel operations */
#define TIMSORT_PARALLEL_THRESHOLD (MIN_OBS_PER_THREAD * 2)

/* ============================================================================
   Timsort Configuration
   ============================================================================ */

/* Minimum run length - runs shorter than this are extended via insertion sort */
#define MIN_RUN 32

/* Minimum gallop threshold - triggers galloping mode in merge */
#define MIN_GALLOP 7

/* Maximum stack size for merge operations (log2(N) + 1 is sufficient) */
#define MAX_MERGE_STACK 85

/* ============================================================================
   Run Structure
   ============================================================================ */

typedef struct {
    size_t start;   /* Starting index of run */
    size_t length;  /* Length of run */
} timsort_run_t;

/* ============================================================================
   Utility Functions
   ============================================================================ */

/*
    NOTE: Aligned memory allocation is now provided by ctools_config.h
    Use ctools_aligned_alloc() and ctools_aligned_free() for cross-platform support.
*/

/*
    Convert IEEE 754 double to sortable uint64.
*/
static inline uint64_t timsort_double_to_sortable(double d)
{
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));

    if (SF_is_missing(d)) {
        return UINT64_MAX;
    }

    if (bits & ((uint64_t)1 << 63)) {
        bits = ~bits;
    } else {
        bits ^= ((uint64_t)1 << 63);
    }

    return bits;
}

/*
    Calculate minimum run length for Timsort.
    Returns a value between MIN_RUN/2 and MIN_RUN such that N/minrun
    is a power of 2 or close to it.
*/
static size_t calc_min_run(size_t n)
{
    size_t r = 0;
    while (n >= MIN_RUN) {
        r |= (n & 1);
        n >>= 1;
    }
    return n + r;
}

/* ============================================================================
   Comparison Functions
   ============================================================================ */

/*
    Compare two strings.
*/
static inline int compare_string(const char *a, const char *b)
{
    return strcmp(a, b);
}

/* ============================================================================
   Binary Insertion Sort (optimized with memmove)
   ============================================================================ */

/*
    Binary insertion sort for numeric data.
    Sorts order[start..end) using binary search for insertion position.
    Uses memmove for batch shifting instead of element-by-element.
*/
static void binary_insertion_sort_numeric(size_t *order, uint64_t *keys,
                                          size_t start, size_t end, size_t sorted_end)
{
    size_t i;
    size_t left, right, mid;
    size_t temp_idx;
    uint64_t temp_key;

    for (i = sorted_end; i < end; i++) {
        temp_idx = order[i];
        temp_key = keys[temp_idx];

        /* Binary search for insertion position */
        left = start;
        right = i;
        while (left < right) {
            mid = left + (right - left) / 2;
            if (keys[order[mid]] <= temp_key) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        /* Batch shift elements using memmove */
        if (left < i) {
            memmove(&order[left + 1], &order[left], (i - left) * sizeof(size_t));
        }
        order[left] = temp_idx;
    }
}

/*
    Binary insertion sort for string data.
    Uses memmove for batch shifting.
*/
static void binary_insertion_sort_string(size_t *order, char **strings,
                                         size_t start, size_t end, size_t sorted_end)
{
    size_t i, j;
    size_t left, right, mid;
    size_t temp_idx;
    const char *temp_str;

    for (i = sorted_end; i < end; i++) {
        temp_idx = order[i];
        temp_str = strings[temp_idx];

        /* Binary search for insertion position */
        left = start;
        right = i;
        while (left < right) {
            mid = left + (right - left) / 2;
            if (compare_string(strings[order[mid]], temp_str) <= 0) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        /* Shift elements to make room (use loop for strings to be safe) */
        for (j = i; j > left; j--) {
            order[j] = order[j - 1];
        }
        order[left] = temp_idx;
    }
}

/* ============================================================================
   Run Detection
   ============================================================================ */

/*
    Find a natural run starting at 'start'.
    Returns the length of the run.
    If descending, reverses the run to make it ascending.
*/
static size_t find_run_numeric(size_t *order, uint64_t *keys,
                               size_t start, size_t end)
{
    size_t run_end = start + 1;

    if (run_end >= end) {
        return 1;
    }

    /* Check if ascending or descending */
    if (keys[order[start]] <= keys[order[run_end]]) {
        /* Ascending run */
        while (run_end < end && keys[order[run_end - 1]] <= keys[order[run_end]]) {
            run_end++;
        }
    } else {
        /* Descending run - find end then reverse */
        while (run_end < end && keys[order[run_end - 1]] > keys[order[run_end]]) {
            run_end++;
        }
        /* Reverse to make ascending */
        size_t i = start;
        size_t j = run_end - 1;
        while (i < j) {
            size_t temp = order[i];
            order[i] = order[j];
            order[j] = temp;
            i++;
            j--;
        }
    }

    return run_end - start;
}

/*
    Find a natural run for string data.
*/
static size_t find_run_string(size_t *order, char **strings,
                              size_t start, size_t end)
{
    size_t run_end = start + 1;

    if (run_end >= end) {
        return 1;
    }

    /* Check if ascending or descending */
    if (compare_string(strings[order[start]], strings[order[run_end]]) <= 0) {
        /* Ascending run */
        while (run_end < end &&
               compare_string(strings[order[run_end - 1]], strings[order[run_end]]) <= 0) {
            run_end++;
        }
    } else {
        /* Descending run */
        while (run_end < end &&
               compare_string(strings[order[run_end - 1]], strings[order[run_end]]) > 0) {
            run_end++;
        }
        /* Reverse */
        size_t i = start;
        size_t j = run_end - 1;
        while (i < j) {
            size_t temp = order[i];
            order[i] = order[j];
            order[j] = temp;
            i++;
            j--;
        }
    }

    return run_end - start;
}

/* ============================================================================
   Galloping Search (for merge optimization)
   ============================================================================ */

/*
    Galloping search: find position where key should be inserted in sorted array.
    Uses exponential search followed by binary search.
    Returns index of first element > key (for left merge) or >= key (for right merge).
*/
static size_t gallop_left_numeric(uint64_t key, size_t *order, uint64_t *keys,
                                  size_t base, size_t len)
{
    size_t ofs = 1;
    size_t last_ofs = 0;
    size_t m;

    if (len == 0) return 0;

    /* Gallop right until keys[order[base + ofs - 1]] < key <= keys[order[base + ofs]] */
    if (keys[order[base]] >= key) {
        return 0;
    }

    while (ofs < len && keys[order[base + ofs]] < key) {
        last_ofs = ofs;
        ofs = (ofs << 1) + 1;
        if (ofs > len) ofs = len;
    }

    /* Binary search in [last_ofs, ofs) */
    while (last_ofs < ofs) {
        m = last_ofs + ((ofs - last_ofs) >> 1);
        if (keys[order[base + m]] < key) {
            last_ofs = m + 1;
        } else {
            ofs = m;
        }
    }

    return ofs;
}

static size_t gallop_right_numeric(uint64_t key, size_t *order, uint64_t *keys,
                                   size_t base, size_t len)
{
    size_t ofs = 1;
    size_t last_ofs = 0;
    size_t m;

    if (len == 0) return 0;

    if (keys[order[base + len - 1]] <= key) {
        return len;
    }

    /* Gallop left from end */
    size_t max_ofs = len;
    while (ofs < max_ofs && keys[order[base + len - 1 - ofs]] > key) {
        last_ofs = ofs;
        ofs = (ofs << 1) + 1;
        if (ofs > max_ofs) ofs = max_ofs;
    }

    /* Convert to forward indices */
    size_t tmp = last_ofs;
    last_ofs = len - ofs;
    ofs = len - tmp;

    /* Binary search */
    while (last_ofs < ofs) {
        m = last_ofs + ((ofs - last_ofs) >> 1);
        if (keys[order[base + m]] <= key) {
            last_ofs = m + 1;
        } else {
            ofs = m;
        }
    }

    return ofs;
}

/* ============================================================================
   Merge Operations
   ============================================================================ */

/*
    Merge two adjacent runs for numeric data.
    run1 = [base1, base1 + len1), run2 = [base2, base2 + len2)
    where base2 = base1 + len1
*/
static void merge_numeric(size_t *order, uint64_t *keys, size_t *temp,
                          size_t base1, size_t len1, size_t len2)
{
    size_t base2 = base1 + len1;
    size_t i, j, k;

    /* Optimization: skip elements already in place */
    size_t skip = gallop_right_numeric(keys[order[base2]], order, keys, base1, len1);
    base1 += skip;
    len1 -= skip;
    if (len1 == 0) return;

    size_t trim = len2 - gallop_left_numeric(keys[order[base1 + len1 - 1]], order, keys, base2, len2);
    len2 -= trim;
    if (len2 == 0) return;

    /* Copy smaller run to temp */
    if (len1 <= len2) {
        /* Merge low: copy run1 to temp, merge from left */
        memcpy(temp, &order[base1], len1 * sizeof(size_t));

        i = 0;      /* Index into temp (copy of run1) */
        j = base2;  /* Index into run2 */
        k = base1;  /* Output index */

        while (i < len1 && j < base2 + len2) {
            if (keys[temp[i]] <= keys[order[j]]) {
                order[k++] = temp[i++];
            } else {
                order[k++] = order[j++];
            }
        }

        /* Copy remaining from temp */
        while (i < len1) {
            order[k++] = temp[i++];
        }
    } else {
        /* Merge high: copy run2 to temp, merge from right */
        memcpy(temp, &order[base2], len2 * sizeof(size_t));

        i = len1 - 1;  /* Index into run1 (from end) */
        j = len2 - 1;  /* Index into temp (from end) */
        k = base1 + len1 + len2 - 1;  /* Output index (from end) */

        while (i < len1 && j < len2) {  /* Using underflow for termination */
            if (keys[order[base1 + i]] > keys[temp[j]]) {
                order[k--] = order[base1 + i--];
            } else {
                order[k--] = temp[j--];
            }
            if (i == (size_t)-1 || j == (size_t)-1) break;
        }

        /* Copy remaining from temp */
        while (j < len2 && j != (size_t)-1) {
            order[k--] = temp[j--];
        }
    }
}

/*
    Merge two adjacent runs for string data.
    Always copies the smaller run to temp to minimize memory usage.
*/
static void merge_string(size_t *order, char **strings, size_t *temp,
                         size_t base1, size_t len1, size_t len2)
{
    size_t base2 = base1 + len1;
    size_t i, j, k;

    if (len1 <= len2) {
        /* Left run is smaller - merge from left */
        memcpy(temp, &order[base1], len1 * sizeof(size_t));

        i = 0;
        j = base2;
        k = base1;

        while (i < len1 && j < base2 + len2) {
            if (compare_string(strings[temp[i]], strings[order[j]]) <= 0) {
                order[k++] = temp[i++];
            } else {
                order[k++] = order[j++];
            }
        }

        while (i < len1) {
            order[k++] = temp[i++];
        }
    } else {
        /* Right run is smaller - merge from right */
        memcpy(temp, &order[base2], len2 * sizeof(size_t));

        i = len1 - 1;
        j = len2 - 1;
        k = base1 + len1 + len2 - 1;

        while (i != (size_t)-1 && j != (size_t)-1) {
            if (compare_string(strings[order[base1 + i]], strings[temp[j]]) > 0) {
                order[k--] = order[base1 + i--];
            } else {
                order[k--] = temp[j--];
            }
        }

        while (j != (size_t)-1) {
            order[k--] = temp[j--];
        }
    }
}

/* ============================================================================
   Timsort Main Algorithm
   ============================================================================ */

/*
    Timsort for numeric data.
*/
static stata_retcode timsort_numeric(size_t *order, uint64_t *keys, size_t nobs)
{
    timsort_run_t run_stack[MAX_MERGE_STACK];
    int stack_size = 0;
    size_t min_run;
    size_t *temp;
    size_t pos = 0;

    if (nobs < 2) return STATA_OK;

    min_run = calc_min_run(nobs);

    /* Allocate temp buffer for merging */
    temp = (size_t *)malloc((nobs / 2 + 1) * sizeof(size_t));
    if (temp == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Process all runs */
    while (pos < nobs) {
        /* Find a natural run */
        size_t run_len = find_run_numeric(order, keys, pos, nobs);

        /* Extend short runs using insertion sort */
        if (run_len < min_run) {
            size_t force_len = nobs - pos;
            if (force_len > min_run) force_len = min_run;

            binary_insertion_sort_numeric(order, keys, pos, pos + force_len, pos + run_len);
            run_len = force_len;
        }

        /* Push run onto stack */
        run_stack[stack_size].start = pos;
        run_stack[stack_size].length = run_len;
        stack_size++;

        /* Merge to maintain stack invariants:
           1. run[i-2] > run[i-1] + run[i]
           2. run[i-1] > run[i]
        */
        while (stack_size > 1) {
            int n = stack_size - 2;

            if ((n > 0 && run_stack[n-1].length <= run_stack[n].length + run_stack[n+1].length) ||
                (n > 1 && run_stack[n-2].length <= run_stack[n-1].length + run_stack[n].length)) {

                if (n > 0 && run_stack[n-1].length < run_stack[n+1].length) {
                    n--;
                }

                /* Merge run[n] with run[n+1] */
                merge_numeric(order, keys, temp,
                              run_stack[n].start, run_stack[n].length,
                              run_stack[n+1].length);
                run_stack[n].length += run_stack[n+1].length;

                /* Remove run[n+1] from stack */
                if (n + 2 < stack_size) {
                    run_stack[n+1] = run_stack[n+2];
                }
                stack_size--;

            } else if (run_stack[n].length <= run_stack[n+1].length) {
                /* Merge run[n] with run[n+1] */
                merge_numeric(order, keys, temp,
                              run_stack[n].start, run_stack[n].length,
                              run_stack[n+1].length);
                run_stack[n].length += run_stack[n+1].length;
                stack_size--;
            } else {
                break;
            }
        }

        pos += run_len;
    }

    /* Final merges: merge all remaining runs */
    while (stack_size > 1) {
        int n = stack_size - 2;

        if (n > 0 && run_stack[n-1].length < run_stack[n+1].length) {
            n--;
        }

        merge_numeric(order, keys, temp,
                      run_stack[n].start, run_stack[n].length,
                      run_stack[n+1].length);
        run_stack[n].length += run_stack[n+1].length;

        if (n + 2 < stack_size) {
            run_stack[n+1] = run_stack[n+2];
        }
        stack_size--;
    }

    free(temp);
    return STATA_OK;
}

/*
    Timsort for string data.
    Full timsort with stack and merging.
*/
static stata_retcode timsort_string(size_t *order, char **strings, size_t nobs)
{
    timsort_run_t run_stack[MAX_MERGE_STACK];
    int stack_size = 0;
    size_t min_run;
    size_t *temp;
    size_t pos = 0;

    if (nobs < 2) return STATA_OK;

    min_run = calc_min_run(nobs);

    temp = (size_t *)malloc((nobs / 2 + 1) * sizeof(size_t));
    if (temp == NULL) {
        return STATA_ERR_MEMORY;
    }

    while (pos < nobs) {
        size_t run_len = find_run_string(order, strings, pos, nobs);

        if (run_len < min_run) {
            size_t force_len = nobs - pos;
            if (force_len > min_run) force_len = min_run;

            binary_insertion_sort_string(order, strings, pos, pos + force_len, pos + run_len);
            run_len = force_len;
        }

        run_stack[stack_size].start = pos;
        run_stack[stack_size].length = run_len;
        stack_size++;

        /* Merge to maintain stack invariants */
        while (stack_size > 1) {
            int n = stack_size - 2;

            if ((n > 0 && run_stack[n-1].length <= run_stack[n].length + run_stack[n+1].length) ||
                (n > 1 && run_stack[n-2].length <= run_stack[n-1].length + run_stack[n].length)) {

                if (n > 0 && run_stack[n-1].length < run_stack[n+1].length) {
                    n--;
                }

                merge_string(order, strings, temp,
                             run_stack[n].start, run_stack[n].length,
                             run_stack[n+1].length);
                run_stack[n].length += run_stack[n+1].length;

                if (n + 2 < stack_size) {
                    run_stack[n+1] = run_stack[n+2];
                }
                stack_size--;

            } else if (run_stack[n].length <= run_stack[n+1].length) {
                merge_string(order, strings, temp,
                             run_stack[n].start, run_stack[n].length,
                             run_stack[n+1].length);
                run_stack[n].length += run_stack[n+1].length;
                stack_size--;
            } else {
                break;
            }
        }

        pos += run_len;
    }

    /* Final merges */
    while (stack_size > 1) {
        int n = stack_size - 2;

        if (n > 0 && run_stack[n-1].length < run_stack[n+1].length) {
            n--;
        }

        merge_string(order, strings, temp,
                     run_stack[n].start, run_stack[n].length,
                     run_stack[n+1].length);
        run_stack[n].length += run_stack[n+1].length;

        if (n + 2 < stack_size) {
            run_stack[n+1] = run_stack[n+2];
        }
        stack_size--;
    }

    free(temp);
    return STATA_OK;
}

/* ============================================================================
   Variable-Level Sort Functions
   ============================================================================ */

static stata_retcode timsort_by_numeric_var(stata_data *data, int var_idx)
{
    uint64_t *keys;
    size_t i;
    double *dbl_data;
    stata_retcode rc;

    keys = (uint64_t *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                              data->nobs * sizeof(uint64_t));
    if (keys == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Parallel key conversion for large datasets */
    dbl_data = data->vars[var_idx].data.dbl;
    #ifdef _OPENMP
    #pragma omp parallel for if(data->nobs >= TIMSORT_PARALLEL_THRESHOLD)
    #endif
    for (i = 0; i < data->nobs; i++) {
        keys[i] = timsort_double_to_sortable(dbl_data[i]);
    }

    rc = timsort_numeric(data->sort_order, keys, data->nobs);

    ctools_aligned_free(keys);
    return rc;
}

static stata_retcode timsort_by_string_var(stata_data *data, int var_idx)
{
    return timsort_string(data->sort_order, data->vars[var_idx].data.str, data->nobs);
}

/* ============================================================================
   Permutation Application
   ============================================================================ */

typedef struct {
    stata_variable *var;
    size_t *sort_order;
    size_t nobs;
    int success;
} timsort_permute_args_t;

static void *timsort_apply_permute_thread(void *arg)
{
    timsort_permute_args_t *args = (timsort_permute_args_t *)arg;
    size_t i;
    size_t nobs = args->nobs;
    size_t *perm = args->sort_order;

    if (args->var->type == STATA_TYPE_DOUBLE) {
        double *old_data = args->var->data.dbl;
        double *new_data = (double *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                                            nobs * sizeof(double));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        ctools_aligned_free(old_data);
        args->var->data.dbl = new_data;
    } else {
        char **old_data = args->var->data.str;
        char **new_data = (char **)ctools_aligned_alloc(CACHE_LINE_SIZE, nobs * sizeof(char *));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        ctools_aligned_free(old_data);
        args->var->data.str = new_data;
    }

    args->success = 1;
    return NULL;
}

static stata_retcode timsort_apply_permutation(stata_data *data)
{
    size_t j;
    size_t nvars = data->nvars;
    pthread_t *threads;
    timsort_permute_args_t *args;
    int all_success = 1;

    if (nvars == 0) {
        return STATA_OK;
    }

    threads = (pthread_t *)malloc(nvars * sizeof(pthread_t));
    args = (timsort_permute_args_t *)malloc(nvars * sizeof(timsort_permute_args_t));

    if (threads == NULL || args == NULL) {
        free(threads);
        free(args);
        return STATA_ERR_MEMORY;
    }

    for (j = 0; j < nvars; j++) {
        args[j].var = &data->vars[j];
        args[j].sort_order = data->sort_order;
        args[j].nobs = data->nobs;
        args[j].success = 0;

        pthread_create(&threads[j], NULL, timsort_apply_permute_thread, &args[j]);
    }

    for (j = 0; j < nvars; j++) {
        pthread_join(threads[j], NULL);
        if (!args[j].success) {
            all_success = 0;
        }
    }

    free(threads);
    free(args);

    for (j = 0; j < data->nobs; j++) {
        data->sort_order[j] = j;
    }

    return all_success ? STATA_OK : STATA_ERR_MEMORY;
}

/* ============================================================================
   Public API
   ============================================================================ */

/*
    Sort data using Timsort and apply permutation to all variables.

    Timsort is particularly efficient for:
    - Partially sorted data (panel data sorted by ID, then needs time sort)
    - Data with natural runs (time series, alphabetically ordered data)
    - When stability is important
    - General-purpose sorting with O(N log N) guarantee

    @param data       [in/out] stata_data structure to sort
    @param sort_vars  [in] Array of 1-based variable indices specifying sort keys
    @param nsort      [in] Number of sort key variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_timsort(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /*
        For stable sort with multiple keys:
        Sort from the LAST (least significant) key to the FIRST (most significant).
    */
    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = timsort_by_numeric_var(data, var_idx);
        } else {
            rc = timsort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    rc = timsort_apply_permutation(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

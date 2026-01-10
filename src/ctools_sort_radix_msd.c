/*
    ctools_sort_radix_msd.c
    Parallel MSD (Most Significant Digit) Radix Sort Module

    MSD radix sort is optimized for variable-length strings where it can
    short-circuit on unique prefixes. Unlike LSD which must process all
    characters, MSD stops as soon as elements are distinguishable.

    Key advantages over LSD for strings:
    - Early termination: stops when all elements in a bucket are identical
    - Better cache locality: strings with common prefixes stay together
    - O(N*k) where k is the average distinguishing prefix length (not max length)

    Algorithm:
    - Process characters from most significant (leftmost) to least
    - Recursively sort each non-empty bucket
    - Use insertion sort for small buckets (better cache behavior)

    Parallelization Strategy:
    - Top-level histogram and scatter are parallelized
    - Recursive calls on independent buckets run in parallel via OpenMP tasks
    - Falls back to sequential insertion sort for small subarrays

    Numeric Sorting:
    - Uses same IEEE 754 bit manipulation as LSD
    - Processes bytes from MSB (byte 7) to LSB (byte 0)
    - Single pass when data has low entropy in high bytes
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/* Threshold for switching to insertion sort */
#define MSD_INSERTION_THRESHOLD 64

/* Minimum bucket size for parallel recursion */
#define MSD_PARALLEL_THRESHOLD 10000

/* ============================================================================
   Utility Functions
   ============================================================================ */

/*
    Aligned memory allocation wrapper.
*/
static void *msd_aligned_alloc(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(__APPLE__) || defined(__linux__)
    if (posix_memalign(&ptr, alignment, size) != 0) {
        return NULL;
    }
    return ptr;
#else
    return malloc(size);
#endif
}

/*
    Convert IEEE 754 double to sortable uint64.
    Same transformation as LSD for consistency.
*/
static inline uint64_t msd_double_to_sortable(double d)
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

/* ============================================================================
   Insertion Sort (for small subarrays)
   ============================================================================ */

/*
    Insertion sort for numeric data using precomputed keys.
    Used when subarray size falls below MSD_INSERTION_THRESHOLD.
*/
static void insertion_sort_numeric(size_t *order, uint64_t *keys,
                                   size_t start, size_t end)
{
    size_t i, j;
    size_t temp_idx;
    uint64_t temp_key;

    for (i = start + 1; i < end; i++) {
        temp_idx = order[i];
        temp_key = keys[temp_idx];
        j = i;

        while (j > start && keys[order[j - 1]] > temp_key) {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp_idx;
    }
}

/*
    Insertion sort for string data.
    Compares strings starting from a given character position.
*/
static void insertion_sort_string(size_t *order, char **strings,
                                  size_t *str_lengths, size_t start,
                                  size_t end, size_t char_pos)
{
    size_t i, j;
    size_t temp_idx;

    for (i = start + 1; i < end; i++) {
        temp_idx = order[i];
        j = i;

        while (j > start) {
            size_t prev_idx = order[j - 1];
            /* Compare strings from char_pos onwards */
            const char *s1 = strings[prev_idx] + char_pos;
            const char *s2 = strings[temp_idx] + char_pos;
            size_t len1 = str_lengths[prev_idx] > char_pos ?
                          str_lengths[prev_idx] - char_pos : 0;
            size_t len2 = str_lengths[temp_idx] > char_pos ?
                          str_lengths[temp_idx] - char_pos : 0;

            /* Compare using available length */
            int cmp;
            size_t min_len = len1 < len2 ? len1 : len2;
            if (min_len == 0) {
                cmp = (int)len1 - (int)len2;
            } else {
                cmp = memcmp(s1, s2, min_len);
                if (cmp == 0) {
                    cmp = (int)len1 - (int)len2;
                }
            }

            if (cmp <= 0) break;

            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp_idx;
    }
}

/* ============================================================================
   MSD Radix Sort for Numeric Data
   ============================================================================ */

/*
    Recursive MSD radix sort for numeric data.
    Processes one byte at a time from MSB to LSB.
*/
static void msd_sort_numeric_recursive(size_t *order, size_t *temp_order,
                                       uint64_t *keys, size_t start,
                                       size_t end, int byte_pos)
{
    size_t counts[RADIX_SIZE] = {0};
    size_t offsets[RADIX_SIZE];
    size_t bucket_starts[RADIX_SIZE];
    size_t i, n;
    uint8_t byte_val;
    int shift;
    int non_empty_buckets = 0;
    int b;

    n = end - start;

    /* Base case: use insertion sort for small arrays */
    if (n <= MSD_INSERTION_THRESHOLD) {
        insertion_sort_numeric(order, keys, start, end);
        return;
    }

    /* Base case: processed all bytes */
    if (byte_pos < 0) {
        return;
    }

    shift = byte_pos * RADIX_BITS;

    /* Count occurrences of each byte value */
    for (i = start; i < end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        counts[byte_val]++;
    }

    /* Check for uniform distribution (all in one bucket) */
    for (b = 0; b < RADIX_SIZE; b++) {
        if (counts[b] > 0) {
            non_empty_buckets++;
            if (non_empty_buckets > 1) break;
        }
    }

    /* If all elements in same bucket, skip this byte and continue */
    if (non_empty_buckets <= 1) {
        msd_sort_numeric_recursive(order, temp_order, keys, start, end,
                                   byte_pos - 1);
        return;
    }

    /* Compute prefix sums (bucket starting positions) */
    offsets[0] = start;
    bucket_starts[0] = start;
    for (b = 1; b < RADIX_SIZE; b++) {
        offsets[b] = offsets[b - 1] + counts[b - 1];
        bucket_starts[b] = offsets[b];
    }

    /* Scatter elements to their bucket positions */
    for (i = start; i < end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        temp_order[offsets[byte_val]++] = order[i];
    }

    /* Copy back to order array */
    memcpy(&order[start], &temp_order[start], n * sizeof(size_t));

    /* Recursively sort each non-empty bucket */
    /* Use OpenMP tasks for large buckets */
    #ifdef _OPENMP
    #pragma omp parallel if(n >= MSD_PARALLEL_THRESHOLD)
    {
        #pragma omp single nowait
        {
    #endif
            for (b = 0; b < RADIX_SIZE; b++) {
                if (counts[b] > 1) {
                    size_t bucket_end = bucket_starts[b] + counts[b];
    #ifdef _OPENMP
                    if (counts[b] >= MSD_PARALLEL_THRESHOLD) {
                        #pragma omp task
                        msd_sort_numeric_recursive(order, temp_order, keys,
                                                   bucket_starts[b], bucket_end,
                                                   byte_pos - 1);
                    } else {
                        msd_sort_numeric_recursive(order, temp_order, keys,
                                                   bucket_starts[b], bucket_end,
                                                   byte_pos - 1);
                    }
    #else
                    msd_sort_numeric_recursive(order, temp_order, keys,
                                               bucket_starts[b], bucket_end,
                                               byte_pos - 1);
    #endif
                }
            }
    #ifdef _OPENMP
        }
        #pragma omp taskwait
    }
    #endif
}

/*
    Sort by a single numeric variable using MSD radix sort.
*/
static stata_retcode msd_sort_by_numeric_var(stata_data *data, int var_idx)
{
    size_t *temp_order;
    uint64_t *keys;
    size_t i;
    double *dbl_data;

    /* Allocate temporary arrays */
    temp_order = (size_t *)msd_aligned_alloc(CACHE_LINE_SIZE,
                                              data->nobs * sizeof(size_t));
    keys = (uint64_t *)msd_aligned_alloc(CACHE_LINE_SIZE,
                                          data->nobs * sizeof(uint64_t));

    if (temp_order == NULL || keys == NULL) {
        free(temp_order);
        free(keys);
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable uint64 keys */
    dbl_data = data->vars[var_idx].data.dbl;
    for (i = 0; i < data->nobs; i++) {
        keys[i] = msd_double_to_sortable(dbl_data[i]);
    }

    /* Sort starting from MSB (byte 7) */
    msd_sort_numeric_recursive(data->sort_order, temp_order, keys,
                               0, data->nobs, 7);

    free(temp_order);
    free(keys);
    return STATA_OK;
}

/* ============================================================================
   MSD Radix Sort for String Data
   ============================================================================ */

/*
    Recursive MSD radix sort for string data.
    Processes one character at a time from leftmost to rightmost.
*/
static void msd_sort_string_recursive(size_t *order, size_t *temp_order,
                                      char **strings, size_t *str_lengths,
                                      size_t start, size_t end, size_t char_pos)
{
    /* 257 buckets: 0 for end-of-string, 1-256 for character values */
    size_t counts[257] = {0};
    size_t offsets[257];
    size_t bucket_starts[257];
    size_t i, n, idx, len;
    unsigned int bucket;
    int non_empty_buckets = 0;
    int b;

    n = end - start;

    /* Base case: use insertion sort for small arrays */
    if (n <= MSD_INSERTION_THRESHOLD) {
        insertion_sort_string(order, strings, str_lengths, start, end, char_pos);
        return;
    }

    /* Count occurrences of each character at current position */
    /* Bucket 0 = end-of-string, buckets 1-256 = character value + 1 */
    for (i = start; i < end; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos >= len) {
            bucket = 0;  /* End of string */
        } else {
            bucket = (unsigned char)strings[idx][char_pos] + 1;
        }
        counts[bucket]++;
    }

    /* Check for uniform distribution */
    for (b = 0; b < 257; b++) {
        if (counts[b] > 0) {
            non_empty_buckets++;
            if (non_empty_buckets > 1) break;
        }
    }

    /* If all elements in same bucket */
    if (non_empty_buckets <= 1) {
        /* If all strings ended, we're done */
        if (counts[0] == n) {
            return;
        }
        /* Otherwise, continue to next character */
        msd_sort_string_recursive(order, temp_order, strings, str_lengths,
                                  start, end, char_pos + 1);
        return;
    }

    /* Compute prefix sums */
    offsets[0] = start;
    bucket_starts[0] = start;
    for (b = 1; b < 257; b++) {
        offsets[b] = offsets[b - 1] + counts[b - 1];
        bucket_starts[b] = offsets[b];
    }

    /* Scatter elements to their bucket positions */
    for (i = start; i < end; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos >= len) {
            bucket = 0;
        } else {
            bucket = (unsigned char)strings[idx][char_pos] + 1;
        }
        temp_order[offsets[bucket]++] = idx;
    }

    /* Copy back to order array */
    memcpy(&order[start], &temp_order[start], n * sizeof(size_t));

    /* Recursively sort each non-empty bucket (skip bucket 0 - strings that ended) */
    #ifdef _OPENMP
    #pragma omp parallel if(n >= MSD_PARALLEL_THRESHOLD)
    {
        #pragma omp single nowait
        {
    #endif
            for (b = 1; b < 257; b++) {  /* Skip bucket 0 */
                if (counts[b] > 1) {
                    size_t bucket_end = bucket_starts[b] + counts[b];
    #ifdef _OPENMP
                    if (counts[b] >= MSD_PARALLEL_THRESHOLD) {
                        #pragma omp task
                        msd_sort_string_recursive(order, temp_order, strings,
                                                  str_lengths, bucket_starts[b],
                                                  bucket_end, char_pos + 1);
                    } else {
                        msd_sort_string_recursive(order, temp_order, strings,
                                                  str_lengths, bucket_starts[b],
                                                  bucket_end, char_pos + 1);
                    }
    #else
                    msd_sort_string_recursive(order, temp_order, strings,
                                              str_lengths, bucket_starts[b],
                                              bucket_end, char_pos + 1);
    #endif
                }
            }
    #ifdef _OPENMP
        }
        #pragma omp taskwait
    }
    #endif
}

/*
    Sort by a single string variable using MSD radix sort.
*/
static stata_retcode msd_sort_by_string_var(stata_data *data, int var_idx)
{
    size_t *temp_order;
    size_t *str_lengths;
    size_t i;
    char **str_data;

    str_data = data->vars[var_idx].data.str;

    /* Pre-cache string lengths */
    str_lengths = (size_t *)malloc(data->nobs * sizeof(size_t));
    if (str_lengths == NULL) {
        return STATA_ERR_MEMORY;
    }

    for (i = 0; i < data->nobs; i++) {
        str_lengths[i] = strlen(str_data[i]);
    }

    /* Allocate temporary order array */
    temp_order = (size_t *)msd_aligned_alloc(CACHE_LINE_SIZE,
                                              data->nobs * sizeof(size_t));
    if (temp_order == NULL) {
        free(str_lengths);
        return STATA_ERR_MEMORY;
    }

    /* Sort starting from position 0 */
    msd_sort_string_recursive(data->sort_order, temp_order, str_data,
                              str_lengths, 0, data->nobs, 0);

    free(temp_order);
    free(str_lengths);
    return STATA_OK;
}

/* ============================================================================
   Permutation Application (shared with LSD)
   ============================================================================ */

/*
    Thread argument structure for parallel permutation application.
*/
typedef struct {
    stata_variable *var;
    size_t *sort_order;
    size_t nobs;
    int success;
} msd_permute_args_t;

/*
    Thread function: Apply permutation to a single variable's data.
*/
static void *msd_apply_permute_thread(void *arg)
{
    msd_permute_args_t *args = (msd_permute_args_t *)arg;
    size_t i;
    size_t nobs = args->nobs;
    size_t *perm = args->sort_order;

    if (args->var->type == STATA_TYPE_DOUBLE) {
        double *old_data = args->var->data.dbl;
        double *new_data = (double *)msd_aligned_alloc(CACHE_LINE_SIZE,
                                                        nobs * sizeof(double));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        free(old_data);
        args->var->data.dbl = new_data;
    } else {
        char **old_data = args->var->data.str;
        char **new_data = (char **)malloc(nobs * sizeof(char *));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        free(old_data);
        args->var->data.str = new_data;
    }

    args->success = 1;
    return NULL;
}

/*
    Apply sort_order permutation to all variables in parallel.
*/
static stata_retcode msd_apply_permutation(stata_data *data)
{
    size_t j;
    size_t nvars = data->nvars;
    pthread_t *threads;
    msd_permute_args_t *args;
    int all_success = 1;

    if (nvars == 0) {
        return STATA_OK;
    }

    threads = (pthread_t *)malloc(nvars * sizeof(pthread_t));
    args = (msd_permute_args_t *)malloc(nvars * sizeof(msd_permute_args_t));

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

        pthread_create(&threads[j], NULL, msd_apply_permute_thread, &args[j]);
    }

    for (j = 0; j < nvars; j++) {
        pthread_join(threads[j], NULL);
        if (!args[j].success) {
            all_success = 0;
        }
    }

    free(threads);
    free(args);

    /* Reset sort_order to identity */
    for (j = 0; j < data->nobs; j++) {
        data->sort_order[j] = j;
    }

    return all_success ? STATA_OK : STATA_ERR_MEMORY;
}

/* ============================================================================
   Public API
   ============================================================================ */

/*
    Sort data using MSD radix sort and apply permutation to all variables.

    MSD radix sort is particularly efficient for:
    - Variable-length strings with common prefixes
    - Data with low entropy in high-order bits/characters
    - Cases where early termination is beneficial

    @param data       [in/out] stata_data structure to sort
    @param sort_vars  [in] Array of 1-based variable indices specifying sort keys
    @param nsort      [in] Number of sort key variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_radix_msd(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /*
        For stable MSD radix sort with multiple keys:
        Sort from the LAST (least significant) key to the FIRST (most significant).
        This ensures that ties in more significant keys are broken by less significant keys.
    */
    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;  /* Convert to 0-based index */

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = msd_sort_by_numeric_var(data, var_idx);
        } else {
            rc = msd_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    /* Apply the permutation to all variables */
    rc = msd_apply_permutation(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

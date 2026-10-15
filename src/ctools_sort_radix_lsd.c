/*
    ctools_sort_radix_lsd.c
    Parallel LSD Radix Sort Module

    Sorts data loaded from Stata using LSD (Least Significant Digit) radix sort.
    This module is designed to be reusable for any sorting operation on data
    already loaded into stata_data structures.

    Algorithm:
    - 8-bit radix (256 buckets per pass), stable
    - Numeric keys and strings of up to CTOOLS_SORT_RECORD_STRING_MAX bytes
      sort together in one record sort (sort_key_columns): every key is read
      as 64-bit words, and (key, row) records move together, so passes read
      and write sequentially
    - Words sort most significant first, and only rows tied on a word are
      sorted by the next one, so a key list stops after its first key when
      that has no ties
    - Within a word, one LSD pass per byte that differs among the keys; many
      records are first split by their most significant differing byte

    Numeric Sorting:
    - IEEE 754 doubles converted to sortable uint64 via bit manipulation
    - Negative numbers: flip all bits (makes them sort before positives)
    - Positive numbers: flip sign bit only
    - Stata missing values: sort to end, in Stata's . < .a < ... < .z order

    String Sorting:
    - Eight bytes per word, big-endian, zero padded: Stata strings hold no
      zero bytes, so shorter strings sort before longer ones they start
    - Key lists with a longer string: one key at a time, last key first, with
      byte-by-byte LSD passes from the rightmost character and pre-cached
      string lengths

    Parallelization:
    - Parallel histogram: each thread counts its chunk, then combine
    - Parallel scatter: each thread writes to pre-computed offsets
    - Threshold: MIN_OBS_PER_THREAD * 2 observations to enable parallel
    - Tied runs of the record sort are sorted one per thread

    Optimizations:
    - Pointer swapping instead of memcpy between passes
    - Early exit for uniform byte distributions (no-op passes)
    - Aligned memory allocation for better cache/SIMD performance
    - Pre-cached string lengths for string sorting

    Post-Sort Permutation:
    - After sorting, applies permutation to all variables in parallel
    - One thread per variable for cache-friendly access pattern
    - Data becomes physically sorted (sequential access for store phase)
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_threads.h"

/*
    NOTE: Aligned memory allocation is now provided by ctools_config.h
    Use ctools_aligned_alloc() and ctools_aligned_free() for cross-platform support.
    On Windows, _aligned_malloc requires _aligned_free - using regular free() causes heap corruption.

    NOTE: Double-to-sortable conversion is now provided by ctools_types.h
    Use ctools_double_to_sortable(d) for sorting doubles.
*/

/* ============================================================================
   Parallel radix sort structures and functions

   Permutation arrays use perm_idx_t (uint32_t) instead of size_t for:
   - 50% memory reduction (4 bytes vs 8 bytes per index)
   - Better cache utilization during sorting
   - Sufficient for Stata's 2^31 observation limit
   ============================================================================ */

/*
    Reusable allocation structure for parallel radix sort.
    Allocated once and reused across all 8 passes to avoid repeated malloc/free.
*/
typedef struct {
    int num_threads;
    size_t **all_local_counts;   /* [num_threads][RADIX_SIZE] */
    size_t *global_counts;       /* [RADIX_SIZE] */
    size_t *global_offsets;      /* [RADIX_SIZE] */
    size_t **thread_offsets;     /* [num_threads][RADIX_SIZE] */
} radix_sort_context_t;

/*
    Allocate reusable context for parallel radix sort.
    Returns NULL on allocation failure.
*/
static radix_sort_context_t *radix_context_alloc(int num_threads)
{
    radix_sort_context_t *ctx;
    int t;

    ctx = (radix_sort_context_t *)calloc(1, sizeof(radix_sort_context_t));
    if (ctx == NULL) return NULL;

    ctx->num_threads = num_threads;
    ctx->all_local_counts = (size_t **)calloc(num_threads, sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
    ctx->thread_offsets = (size_t **)calloc(num_threads, sizeof(size_t *));

    if (!ctx->all_local_counts || !ctx->global_counts ||
        !ctx->global_offsets || !ctx->thread_offsets) {
        goto fail;
    }

    /* Initialize pointer arrays to NULL for safe cleanup */
    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = NULL;
        ctx->thread_offsets[t] = NULL;
    }

    /* Allocate per-thread arrays */
    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
        ctx->thread_offsets[t] = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
        if (!ctx->all_local_counts[t] || !ctx->thread_offsets[t]) {
            goto fail;
        }
    }

    return ctx;

fail:
    if (ctx) {
        if (ctx->all_local_counts) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->all_local_counts[t]);
            }
        }
        if (ctx->thread_offsets) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->thread_offsets[t]);
            }
        }
        free(ctx->all_local_counts);
        free(ctx->global_counts);
        free(ctx->global_offsets);
        free(ctx->thread_offsets);
        free(ctx);
    }
    return NULL;
}

/*
    Free radix sort context.
*/
static void radix_context_free(radix_sort_context_t *ctx)
{
    int t;
    if (ctx == NULL) return;

    for (t = 0; t < ctx->num_threads; t++) {
        free(ctx->all_local_counts[t]);
        free(ctx->thread_offsets[t]);
    }
    free(ctx->all_local_counts);
    free(ctx->global_counts);
    free(ctx->global_offsets);
    free(ctx->thread_offsets);
    free(ctx);
}

/* ============================================================================
   String sorting with pre-cached lengths and parallel histogram/scatter
   ============================================================================ */

/*
    String sort context for reusable allocations.
*/
typedef struct {
    int num_threads;
    size_t **all_local_counts;  /* [num_threads][RADIX_SIZE+1] */
    size_t *global_counts;      /* [RADIX_SIZE+1] */
    size_t *global_offsets;     /* [RADIX_SIZE+1] */
    size_t **thread_offsets;    /* [num_threads][RADIX_SIZE+1] */
} string_sort_context_t;

static string_sort_context_t *string_context_alloc(int num_threads)
{
    string_sort_context_t *ctx;
    int t;

    ctx = (string_sort_context_t *)calloc(1, sizeof(string_sort_context_t));
    if (ctx == NULL) return NULL;

    ctx->num_threads = num_threads;
    ctx->all_local_counts = (size_t **)calloc(num_threads, sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE + 1, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc((RADIX_SIZE + 1) * sizeof(size_t));
    ctx->thread_offsets = (size_t **)calloc(num_threads, sizeof(size_t *));

    if (!ctx->all_local_counts || !ctx->global_counts ||
        !ctx->global_offsets || !ctx->thread_offsets) {
        goto fail;
    }

    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = NULL;
        ctx->thread_offsets[t] = NULL;
    }

    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = (size_t *)calloc(RADIX_SIZE + 1, sizeof(size_t));
        ctx->thread_offsets[t] = (size_t *)malloc((RADIX_SIZE + 1) * sizeof(size_t));
        if (!ctx->all_local_counts[t] || !ctx->thread_offsets[t]) {
            goto fail;
        }
    }

    return ctx;

fail:
    if (ctx) {
        if (ctx->all_local_counts) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->all_local_counts[t]);
            }
        }
        if (ctx->thread_offsets) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->thread_offsets[t]);
            }
        }
        free(ctx->all_local_counts);
        free(ctx->global_counts);
        free(ctx->global_offsets);
        free(ctx->thread_offsets);
        free(ctx);
    }
    return NULL;
}

static void string_context_free(string_sort_context_t *ctx)
{
    int t;
    if (ctx == NULL) return;

    for (t = 0; t < ctx->num_threads; t++) {
        free(ctx->all_local_counts[t]);
        free(ctx->thread_offsets[t]);
    }
    free(ctx->all_local_counts);
    free(ctx->global_counts);
    free(ctx->global_offsets);
    free(ctx->thread_offsets);
    free(ctx);
}

/*
    Parallel string radix sort pass.
    Returns 1 if skipped (uniform), 0 otherwise.
*/
static int radix_sort_pass_string_parallel(perm_idx_t *order,
                                           perm_idx_t *temp_order,
                                           char **strings,
                                           size_t *str_lengths,
                                           size_t nobs,
                                           size_t char_pos,
                                           string_sort_context_t *ctx)
{
    size_t chunk_size;
    int t, b;
    int num_threads = ctx->num_threads;

    chunk_size = (nobs + num_threads - 1) / num_threads;

    /* Phase 1: Parallel histogram (OpenMP) */
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t *my_counts = ctx->all_local_counts[tid];
        memset(my_counts, 0, (RADIX_SIZE + 1) * sizeof(size_t));

        for (size_t i = my_start; i < my_end; i++) {
            perm_idx_t idx = order[i];
            size_t len = str_lengths[idx];
            unsigned char byte_val;
            if (char_pos < len && strings[idx] != NULL) {
                byte_val = (unsigned char)strings[idx][char_pos];
            } else {
                byte_val = 0;
            }
            my_counts[byte_val]++;
        }
    }

    /* Combine histograms */
    /* Optimization: iterate bucket-first for better cache locality */
    memset(ctx->global_counts, 0, (RADIX_SIZE + 1) * sizeof(size_t));
    for (b = 0; b <= RADIX_SIZE; b++) {
        for (t = 0; t < num_threads; t++) {
            ctx->global_counts[b] += ctx->all_local_counts[t][b];
        }
    }

    /* Check for uniform distribution */
    {
        int non_empty = 0;
        for (b = 0; b <= RADIX_SIZE; b++) {
            if (ctx->global_counts[b] > 0) {
                non_empty++;
                if (non_empty > 1) break;
            }
        }
        if (non_empty <= 1) {
            return 1;
        }
    }

    /* Compute prefix sums */
    ctx->global_offsets[0] = 0;
    for (b = 1; b <= RADIX_SIZE; b++) {
        ctx->global_offsets[b] = ctx->global_offsets[b - 1] + ctx->global_counts[b - 1];
    }

    /* Compute per-thread offsets */
    for (b = 0; b <= RADIX_SIZE; b++) {
        size_t offset = ctx->global_offsets[b];
        for (t = 0; t < num_threads; t++) {
            ctx->thread_offsets[t][b] = offset;
            offset += ctx->all_local_counts[t][b];
        }
    }

    /* Phase 2: Parallel scatter (OpenMP) */
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t local_offsets[RADIX_SIZE + 1];
        memcpy(local_offsets, ctx->thread_offsets[tid], (RADIX_SIZE + 1) * sizeof(size_t));

        for (size_t i = my_start; i < my_end; i++) {
            perm_idx_t idx = order[i];
            size_t len = str_lengths[idx];
            unsigned char byte_val;
            if (char_pos < len && strings[idx] != NULL) {
                byte_val = (unsigned char)strings[idx][char_pos];
            } else {
                byte_val = 0;
            }
            temp_order[local_offsets[byte_val]++] = idx;
        }
    }

    return 0;
}

/*
    Sequential string radix sort pass with pre-cached lengths.
    Returns 1 if skipped, 0 otherwise.
*/
static int radix_sort_pass_string(perm_idx_t *order,
                                  perm_idx_t *temp_order,
                                  char **strings,
                                  size_t *str_lengths,
                                  size_t nobs,
                                  size_t char_pos)
{
    size_t counts[RADIX_SIZE + 1] = {0};
    size_t offsets[RADIX_SIZE + 1];
    size_t i, idx, len;
    unsigned char byte_val;
    int non_empty = 0;

    /* Count occurrences using pre-cached lengths */
    for (i = 0; i < nobs; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos < len && strings[idx] != NULL) {
            byte_val = (unsigned char)strings[idx][char_pos];
        } else {
            byte_val = 0;
        }
        counts[byte_val]++;
    }

    /* Check for uniform distribution */
    for (i = 0; i <= RADIX_SIZE; i++) {
        if (counts[i] > 0) {
            non_empty++;
            if (non_empty > 1) break;
        }
    }
    if (non_empty <= 1) {
        return 1;
    }

    /* Compute prefix sums */
    offsets[0] = 0;
    for (i = 1; i <= RADIX_SIZE; i++) {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }

    /* Scatter */
    for (i = 0; i < nobs; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos < len && strings[idx] != NULL) {
            byte_val = (unsigned char)strings[idx][char_pos];
        } else {
            byte_val = 0;
        }
        temp_order[offsets[byte_val]++] = idx;
    }

    return 0;
}

/* ============================================================================
   Numeric and string keys: stable radix sort of (key, row) records

   Every key variable reads as a sequence of 64-bit words that compare as
   the variable does: a double is one word (ctools_double_to_sortable), a
   string is its bytes eight at a time, big-endian and zero padded. Stata
   strings hold no zero bytes, so the padding puts a string before every
   longer string that starts with it, and unsigned word order is Stata's
   bytewise string order.

   Records sort one word at a time, most significant first: every row by the
   first word, then each run of rows that tie on it by the next word, and so
   on. A run without ties is finished and later words are never read, so a
   key list stops after its first variable when that has no ties. Within a
   word, (key, row) records move together through stable radix passes that
   read and write memory sequentially, one per byte on which the keys
   differ. Many records are first split by their most significant
   differing byte, so that each bucket sorts only on the bytes that still
   differ within it (often none) and fits in cache. Runs that tie on a word
   are independent, so each sorts on one thread.
   ============================================================================ */

typedef struct {
    uint64_t *keys;     /* sortable key of each record */
    perm_idx_t *rows;   /* zero-based row of each record */
} radix_records_t;

/* The key variables, most significant first. */
typedef struct {
    const stata_variable *vars;
    size_t nvars;
} radix_key_list_t;

/* Runs at most this long sort by insertion instead of radix passes. */
#define RADIX_RUN_INSERTION 24

/* Longer runs split by their most significant differing byte before any
   LSD pass, so that the passes run on buckets that stay in cache. */
#define RADIX_RUN_CACHE 131072

/* Parallel sorts of a word split by its top byte when the keys differ in at
   least this many bytes; with fewer, plain passes move records no more often. */
#define RADIX_SPLIT_MIN_BYTES 3

/* Keys sampled to predict whether such a split pays off. */
#define RADIX_SPLIT_SAMPLE 1024

static int radix_record_threads(size_t nobs)
{
    int threads = ctools_get_openmp_threads();
    if (nobs < (size_t)MIN_OBS_PER_THREAD * 2) return 1;
    if (nobs / MIN_OBS_PER_THREAD < (size_t)threads) threads = (int)(nobs / MIN_OBS_PER_THREAD);
    if (threads > 64) threads = 64;
    return threads < 1 ? 1 : threads;
}

/* Contiguous share of [0, nobs) for worker tid of a team. */
static void radix_record_share(size_t nobs, int team, int tid, size_t *start, size_t *end)
{
    size_t base = nobs / (size_t)team, extra = nobs % (size_t)team;
    *start = (size_t)tid * base + ((size_t)tid < extra ? (size_t)tid : extra);
    *end = *start + base + ((size_t)tid < extra);
}

/*
    Eight bytes of a string from byte `off` on, big-endian, zero padded after
    its end. The string must hold at least `off` bytes; NULL reads as empty.
    Nothing after the terminator is read: once a zero byte is loaded the
    pointer stops advancing, so every later load reads that byte again and
    the loop has no branch that depends on the string's length.
*/
static inline uint64_t radix_string_word(const char *s, size_t off)
{
    static const unsigned char empty = 0;
    const unsigned char *p = s != NULL ? (const unsigned char *)s + off : &empty;
    uint64_t word = 0;
    for (int b = 0; b < 8; b++) {
        const unsigned char c = *p;
        word |= (uint64_t)c << (56 - 8 * b);
        p += c != 0;
    }
    return word;
}

/*
    keys[i] = the word at byte `off` of variable `var` in row rows[i], for
    i < n. Returns the bits on which some key differs from `first`.
*/
static uint64_t radix_gather_words(uint64_t *keys, const perm_idx_t *rows, size_t n,
                                   const stata_variable *var, size_t off, uint64_t first)
{
    uint64_t differ = 0;
    if (var->type == STATA_TYPE_DOUBLE) {
        const double *values = var->data.dbl;
        for (size_t i = 0; i < n; i++) {
            if (i + PREFETCH_DISTANCE < n) CTOOLS_PREFETCH(&values[rows[i + PREFETCH_DISTANCE]]);
            const uint64_t key = ctools_double_to_sortable(values[rows[i]]);
            keys[i] = key;
            differ |= key ^ first;
        }
    } else {
        /* Two loads per row, the pointer and then its bytes: fetch pointers
           twice as far ahead as the bytes they point to. */
        char *const *strings = var->data.str;
        for (size_t i = 0; i < n; i++) {
            if (i + 2 * PREFETCH_DISTANCE < n) {
                CTOOLS_PREFETCH(&strings[rows[i + 2 * PREFETCH_DISTANCE]]);
            }
            if (i + PREFETCH_DISTANCE < n) {
                const char *ahead = strings[rows[i + PREFETCH_DISTANCE]];
                if (ahead != NULL) CTOOLS_PREFETCH(ahead + off);
            }
            const uint64_t key = radix_string_word(strings[rows[i]], off);
            keys[i] = key;
            differ |= key ^ first;
        }
    }
    return differ;
}

/* The first word of rows[0], against which radix_gather_words compares. */
static uint64_t radix_first_word(const perm_idx_t *rows, const stata_variable *var, size_t off)
{
    return var->type == STATA_TYPE_DOUBLE ? ctools_double_to_sortable(var->data.dbl[rows[0]])
                                          : radix_string_word(var->data.str[rows[0]], off);
}

/*
    The word that orders rows tied on `word`, the word at byte *off of
    variable *k: the string's next eight bytes unless it ended within this
    word, else the first word of the next variable. Returns 0 when the tie
    is final.
*/
static int radix_next_word(const radix_key_list_t *list, size_t *k, size_t *off, uint64_t word)
{
    if (list->vars[*k].type == STATA_TYPE_STRING && (word & RADIX_MASK) != 0) {
        *off += 8;
        return 1;
    }
    *off = 0;
    return ++*k < list->nvars;
}

/* 1 when no tie on a word of variable k at byte off can need a later word. */
static int radix_word_is_last(const radix_key_list_t *list, size_t k)
{
    return k + 1 == list->nvars && list->vars[k].type == STATA_TYPE_DOUBLE;
}

/*
    Stable sort of the n records (keys, rows) by key over the bytes set in
    `differ`, with (tkeys, trows) as scratch of the same size: one LSD pass
    per differing byte, after radix_split_run for long runs. Returns 1 when
    the sorted records are in the scratch arrays. Single-threaded.
*/
static int radix_sort_run(uint64_t *keys, perm_idx_t *rows, uint64_t *tkeys, perm_idx_t *trows,
                          size_t n, uint64_t differ, uint32_t (*counts)[RADIX_SIZE]);

/*
    radix_sort_run for a long run with several differing bytes: distribute
    the records into the scratch arrays by the most significant differing
    byte (at `shift`), then sort each bucket over the lower bytes that differ
    within it. The result is in the scratch arrays.
*/
static void radix_split_run(uint64_t *keys, perm_idx_t *rows, uint64_t *tkeys, perm_idx_t *trows,
                            size_t n, int shift, uint32_t (*counts)[RADIX_SIZE])
{
    const uint64_t below = UINT64_MAX >> (64 - shift);
    uint32_t start[RADIX_SIZE + 1] = {0};
    uint64_t bucket_or[RADIX_SIZE], bucket_and[RADIX_SIZE];

    for (int b = 0; b < RADIX_SIZE; b++) {
        bucket_or[b] = 0;
        bucket_and[b] = UINT64_MAX;
    }
    for (size_t i = 0; i < n; i++) {
        const uint64_t key = keys[i];
        const int b = (int)((key >> shift) & RADIX_MASK);
        start[b + 1]++;
        bucket_or[b] |= key;
        bucket_and[b] &= key;
    }
    for (int b = 0; b < RADIX_SIZE; b++) start[b + 1] += start[b];
    {
        uint32_t place[RADIX_SIZE];
        memcpy(place, start, sizeof(place));
        for (size_t i = 0; i < n; i++) {
            const uint64_t key = keys[i];
            const uint32_t dest = place[(key >> shift) & RADIX_MASK]++;
            tkeys[dest] = key;
            trows[dest] = rows[i];
        }
    }
    for (int b = 0; b < RADIX_SIZE; b++) {
        const size_t first = start[b], count = start[b + 1] - first;
        const uint64_t bucket_differ = (bucket_or[b] ^ bucket_and[b]) & below;
        if (count < 2 || bucket_differ == 0) continue;
        if (radix_sort_run(tkeys + first, trows + first, keys + first, rows + first, count,
                           bucket_differ, counts)) {
            memcpy(tkeys + first, keys + first, count * sizeof(uint64_t));
            memcpy(trows + first, rows + first, count * sizeof(perm_idx_t));
        }
    }
}

static int radix_sort_run(uint64_t *keys, perm_idx_t *rows, uint64_t *tkeys, perm_idx_t *trows,
                          size_t n, uint64_t differ, uint32_t (*counts)[RADIX_SIZE])
{
    int shift[8], passes = 0;

    if (differ == 0) return 0;
    if (n <= RADIX_RUN_INSERTION) {
        for (size_t i = 1; i < n; i++) {
            const uint64_t key = keys[i];
            const perm_idx_t row = rows[i];
            size_t j = i;
            for (; j > 0 && keys[j - 1] > key; j--) {
                keys[j] = keys[j - 1];
                rows[j] = rows[j - 1];
            }
            keys[j] = key;
            rows[j] = row;
        }
        return 0;
    }

    for (int b = 0; b < 8; b++) {
        if ((differ >> (b * RADIX_BITS)) & RADIX_MASK) shift[passes++] = b * RADIX_BITS;
    }
    if (n > RADIX_RUN_CACHE && passes > 1) {
        radix_split_run(keys, rows, tkeys, trows, n, shift[passes - 1], counts);
        return 1;
    }
    memset(counts, 0, (size_t)passes * sizeof(*counts));
    for (size_t i = 0; i < n; i++) {
        const uint64_t key = keys[i];
        for (int p = 0; p < passes; p++) counts[p][(key >> shift[p]) & RADIX_MASK]++;
    }
    for (int p = 0; p < passes; p++) {
        uint32_t next = 0, *place = counts[p];
        for (int b = 0; b < RADIX_SIZE; b++) {
            const uint32_t count = place[b];
            place[b] = next;
            next += count;
        }
        for (size_t i = 0; i < n; i++) {
            const uint64_t key = keys[i];
            const uint32_t dest = place[(key >> shift[p]) & RADIX_MASK]++;
            tkeys[dest] = key;
            trows[dest] = rows[i];
        }
        uint64_t *swap_keys = keys;
        perm_idx_t *swap_rows = rows;
        keys = tkeys;
        rows = trows;
        tkeys = swap_keys;
        trows = swap_rows;
    }
    return passes & 1;
}

/*
    Sort n records whose rows are set in cur by the key list from word (k, off)
    on, using tmp as scratch, and write their final order to out, a separate
    array. All three hold the run at their start. Single-threaded.
*/
static void radix_refine_run(const radix_key_list_t *list, size_t k, size_t off,
                             radix_records_t cur, radix_records_t tmp, perm_idx_t *out, size_t n,
                             uint32_t (*counts)[RADIX_SIZE])
{
    for (;;) {
        const stata_variable *var = &list->vars[k];
        const uint64_t differ = radix_gather_words(cur.keys, cur.rows, n, var, off,
                                                   radix_first_word(cur.rows, var, off));
        if (radix_sort_run(cur.keys, cur.rows, tmp.keys, tmp.rows, n, differ, counts)) {
            const radix_records_t swap = cur;
            cur = tmp;
            tmp = swap;
        }
        memcpy(out, cur.rows, n * sizeof(perm_idx_t));
        if (radix_word_is_last(list, k)) return;

        /* One run: continue with the next word without recursing. */
        if (cur.keys[0] == cur.keys[n - 1]) {
            if (!radix_next_word(list, &k, &off, cur.keys[0])) return;
            continue;
        }

        for (size_t i = 0, j; i < n; i = j) {
            const uint64_t key = cur.keys[i];
            for (j = i + 1; j < n && cur.keys[j] == key; j++) {}
            size_t k2 = k, off2 = off;
            if (j - i > 1 && radix_next_word(list, &k2, &off2, key)) {
                const radix_records_t sub = {cur.keys + i, cur.rows + i};
                const radix_records_t sub_tmp = {tmp.keys + i, tmp.rows + i};
                radix_refine_run(list, k2, off2, sub, sub_tmp, out + i, j - i, counts);
            }
        }
        return;
    }
}

/*
    Sort `cur` by key (stable), one pass per byte set in `differ`, using `tmp`
    as scratch of the same size, with at most ctx->num_threads workers.
    Returns the records that hold the result: either `cur` or `tmp`.
*/
static const radix_records_t *radix_sort_records(const radix_records_t *cur,
                                                 const radix_records_t *tmp, size_t nobs,
                                                 uint64_t differ, radix_sort_context_t *ctx)
{
    const int threads = ctx->num_threads;
    (void)threads;  /* only named in OpenMP clauses */

    for (int byte_pos = 0; byte_pos < 8; byte_pos++) {
        const int shift = byte_pos * RADIX_BITS;
        if (((differ >> shift) & RADIX_MASK) == 0) continue;
        const uint64_t *keys = cur->keys;
        const perm_idx_t *rows = cur->rows;
        uint64_t *keys_out = tmp->keys;
        perm_idx_t *rows_out = tmp->rows;

        #pragma omp parallel num_threads(threads) if(threads > 1)
        {
            /* The runtime may grant fewer workers than requested. */
            int team = 1, tid = 0;
            #ifdef _OPENMP
            team = omp_get_num_threads();
            tid = omp_get_thread_num();
            #endif
            size_t start, end;
            radix_record_share(nobs, team, tid, &start, &end);
            size_t *counts = ctx->all_local_counts[tid];
            memset(counts, 0, RADIX_SIZE * sizeof(size_t));
            for (size_t i = start; i < end; i++) counts[(keys[i] >> shift) & RADIX_MASK]++;

            /* Bucket b of worker t starts after every smaller bucket and after
               the same bucket of earlier workers, which keeps the pass stable. */
            #pragma omp barrier
            #pragma omp single
            {
                size_t next = 0;
                for (int b = 0; b < RADIX_SIZE; b++) {
                    for (int t = 0; t < team; t++) {
                        ctx->thread_offsets[t][b] = next;
                        next += ctx->all_local_counts[t][b];
                    }
                }
            }

            size_t *place = ctx->thread_offsets[tid];
            for (size_t i = start; i < end; i++) {
                uint64_t key = keys[i];
                size_t dest = place[(key >> shift) & RADIX_MASK]++;
                keys_out[dest] = key;
                rows_out[dest] = rows[i];
            }
        }
        const radix_records_t *swap = cur;
        cur = tmp;
        tmp = swap;
    }
    return cur;
}

/* The part of records r that starts at record `start`. */
static radix_records_t radix_part(const radix_records_t *r, size_t start)
{
    const radix_records_t part = {r->keys + start, r->rows + start};
    return part;
}

/* dst[0, n) = src[0, n), with every worker of ctx. */
static void radix_copy_records(const radix_records_t *dst, const radix_records_t *src, size_t n,
                               radix_sort_context_t *ctx)
{
    const int threads = ctx->num_threads;
    (void)threads;  /* only named in OpenMP clauses */
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
    for (size_t i = 0; i < n; i++) {
        dst->keys[i] = src->keys[i];
        dst->rows[i] = src->rows[i];
    }
}

/*
    1 when, among `sample` keys spread evenly over keys[0, n), more than
    `limit` fall in one bucket of the byte at `shift` whose keys still differ
    below it: splitting by that byte would leave most records to sort again.
*/
static int radix_split_unproductive(const uint64_t *keys, size_t n, int shift, size_t sample,
                                    size_t limit)
{
    const uint64_t below = UINT64_MAX >> (64 - shift);
    uint64_t key_or[RADIX_SIZE], key_and[RADIX_SIZE];
    size_t count[RADIX_SIZE] = {0};
    const size_t step = n / sample;

    if (step == 0) return 0;
    for (int b = 0; b < RADIX_SIZE; b++) {
        key_or[b] = 0;
        key_and[b] = UINT64_MAX;
    }
    for (size_t i = 0; i < sample; i++) {
        const uint64_t key = keys[i * step];
        const int b = (int)((key >> shift) & RADIX_MASK);
        count[b]++;
        key_or[b] |= key;
        key_and[b] &= key;
    }
    for (int b = 0; b < RADIX_SIZE; b++) {
        if (count[b] > limit && ((key_or[b] ^ key_and[b]) & below) != 0) return 1;
    }
    return 0;
}

/* radix_sort_records, with the result moved to dst unless dst is NULL. */
static const radix_records_t *radix_sort_plain(const radix_records_t *src,
                                               const radix_records_t *scratch,
                                               const radix_records_t *dst, size_t n,
                                               uint64_t differ, radix_sort_context_t *ctx)
{
    const radix_records_t *sorted = radix_sort_records(src, scratch, n, differ, ctx);
    if (dst == NULL || sorted->keys == dst->keys) return sorted;
    radix_copy_records(dst, sorted, n, ctx);
    return dst;
}

/*
    Sort the n records of `src` by key (stable) over the bytes set in `differ`,
    using `scratch` as the other buffer of the same size. The result goes to
    `dst` (src or scratch), or, when dst is NULL, to whichever of them the
    passes end in. Returns the records that hold the result.

    Rather than one pass over every record per differing byte, records are
    first distributed by their most significant differing byte. Each bucket
    then sorts on its own over the lower bytes that differ within it (often
    none or one): buckets of at least `large` records the same way with every
    worker, smaller ones one per thread, in cache, each finishing in dst.
    Keys that differ in few bytes, or whose largest unsorted bucket would
    hold most records, take plain passes instead.
*/
static const radix_records_t *radix_sort_word(const radix_records_t *src,
                                              const radix_records_t *scratch,
                                              const radix_records_t *dst, size_t n,
                                              uint64_t differ, size_t large,
                                              radix_sort_context_t *ctx)
{
    const int threads = ctx->num_threads;
    int top = 7, bytes = 0;
    for (int b = 0; b < 8; b++) bytes += ((differ >> (b * RADIX_BITS)) & RADIX_MASK) != 0;
    if (bytes < RADIX_SPLIT_MIN_BYTES) return radix_sort_plain(src, scratch, dst, n, differ, ctx);
    while (((differ >> (top * RADIX_BITS)) & RADIX_MASK) == 0) top--;
    const int shift = top * RADIX_BITS;
    const uint64_t below = UINT64_MAX >> (64 - shift);
    if (radix_split_unproductive(src->keys, n, shift, RADIX_SPLIT_SAMPLE, RADIX_SPLIT_SAMPLE * 3 / 5)) {
        return radix_sort_plain(src, scratch, dst, n, differ, ctx);
    }

    /* Count each share's buckets, and the bits that differ within each bucket
       as the OR and AND of its keys. Shares are fixed, so the scatter can
       reuse the counts whatever team the runtime grants. */
    uint64_t bucket_or[RADIX_SIZE], bucket_and[RADIX_SIZE];
    size_t bucket_start[RADIX_SIZE + 1];
    for (int b = 0; b < RADIX_SIZE; b++) {
        bucket_or[b] = 0;
        bucket_and[b] = UINT64_MAX;
    }
    const uint64_t *keys = src->keys;
    const perm_idx_t *rows = src->rows;
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
    for (int t = 0; t < threads; t++) {
        size_t start, end;
        radix_record_share(n, threads, t, &start, &end);
        size_t *counts = ctx->all_local_counts[t];
        uint64_t key_or[RADIX_SIZE], key_and[RADIX_SIZE];
        memset(counts, 0, RADIX_SIZE * sizeof(size_t));
        for (int b = 0; b < RADIX_SIZE; b++) {
            key_or[b] = 0;
            key_and[b] = UINT64_MAX;
        }
        for (size_t i = start; i < end; i++) {
            const uint64_t key = keys[i];
            const int b = (int)((key >> shift) & RADIX_MASK);
            counts[b]++;
            key_or[b] |= key;
            key_and[b] &= key;
        }
        #pragma omp critical(radix_sort_word_buckets)
        for (int b = 0; b < RADIX_SIZE; b++) {
            bucket_or[b] |= key_or[b];
            bucket_and[b] &= key_and[b];
        }
    }

    /* Stable placement, as in radix_sort_records. */
    size_t next = 0, largest = 0;
    for (int b = 0; b < RADIX_SIZE; b++) {
        bucket_start[b] = next;
        for (int t = 0; t < threads; t++) {
            ctx->thread_offsets[t][b] = next;
            next += ctx->all_local_counts[t][b];
        }
        if ((bucket_or[b] ^ bucket_and[b]) & below && next - bucket_start[b] > largest) {
            largest = next - bucket_start[b];
        }
    }
    bucket_start[RADIX_SIZE] = next;

    /* Splitting off a minority only adds a pass over the bucket that still
       holds most records and needs sorting (the sample can miss that). */
    if (largest > n / 2) return radix_sort_plain(src, scratch, dst, n, differ, ctx);
    if (dst == NULL) dst = scratch;

    uint64_t *keys_out = scratch->keys;
    perm_idx_t *rows_out = scratch->rows;
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
    for (int t = 0; t < threads; t++) {
        size_t start, end;
        radix_record_share(n, threads, t, &start, &end);
        size_t *place = ctx->thread_offsets[t];
        for (size_t i = start; i < end; i++) {
            const uint64_t key = keys[i];
            const size_t dest = place[(key >> shift) & RADIX_MASK]++;
            keys_out[dest] = key;
            rows_out[dest] = rows[i];
        }
    }

    /* Every bucket is in scratch now; its part of src is free. */
    for (int b = 0; b < RADIX_SIZE; b++) {
        const size_t start = bucket_start[b], count = bucket_start[b + 1] - start;
        const uint64_t bucket_differ = (bucket_or[b] ^ bucket_and[b]) & below;
        if (count < large) continue;
        const radix_records_t from = radix_part(scratch, start), other = radix_part(src, start);
        const radix_records_t to = radix_part(dst, start);
        if (bucket_differ != 0) {
            radix_sort_word(&from, &other, &to, count, bucket_differ, large, ctx);
        } else if (dst != scratch) {
            radix_copy_records(&to, &from, count, ctx);
        }
    }
    #pragma omp parallel num_threads(threads) if(threads > 1)
    {
        uint32_t counts[8][RADIX_SIZE];
        #pragma omp for schedule(dynamic, 1)
        for (int b = 0; b < RADIX_SIZE; b++) {
            const size_t start = bucket_start[b], count = bucket_start[b + 1] - start;
            const uint64_t bucket_differ = (bucket_or[b] ^ bucket_and[b]) & below;
            if (count == 0 || count >= large) continue;
            radix_records_t at = radix_part(scratch, start), other = radix_part(src, start);
            if (bucket_differ != 0 &&
                radix_sort_run(at.keys, at.rows, other.keys, other.rows, count, bucket_differ,
                               counts)) {
                at = other;
            }
            if (at.keys != dst->keys + start) {
                memcpy(dst->keys + start, at.keys, count * sizeof(uint64_t));
                memcpy(dst->rows + start, at.rows, count * sizeof(perm_idx_t));
            }
        }
    }
    return dst;
}

/* A run of tied records to sort by later words. */
typedef struct {
    perm_idx_t start;
    perm_idx_t len;
} radix_run_t;

/*
    Runs of keys[0, n) longer than one record whose tie a later word can
    break: *runs (NULL when there are none) and their number. The scan is
    split into `threads` shares; each owns the runs that start in it.
*/
static stata_retcode radix_collect_runs(const radix_key_list_t *list, size_t k, size_t off,
                                        const uint64_t *keys, size_t n, int threads,
                                        radix_run_t **runs, size_t *nruns)
{
    size_t *found = (size_t *)calloc((size_t)threads + 1, sizeof(size_t));
    radix_run_t *list_out = NULL;
    stata_retcode rc = STATA_OK;

    *runs = NULL;
    *nruns = 0;
    if (!found) return STATA_ERR_MEMORY;

    /* Pass 0 counts the runs of each share, pass 1 writes them after the
       runs of earlier shares. */
    for (int pass = 0; pass < 2; pass++) {
        #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
        for (int t = 0; t < threads; t++) {
            size_t start, end, count = 0;
            radix_record_share(n, threads, t, &start, &end);
            if (start > 0) {
                while (start < end && keys[start] == keys[start - 1]) start++;
            }
            for (size_t i = start, j; i < end; i = j) {
                const uint64_t key = keys[i];
                for (j = i + 1; j < n && keys[j] == key; j++) {}
                size_t k2 = k, off2 = off;
                if (j - i > 1 && radix_next_word(list, &k2, &off2, key)) {
                    if (pass == 1) {
                        list_out[found[t] + count].start = (perm_idx_t)i;
                        list_out[found[t] + count].len = (perm_idx_t)(j - i);
                    }
                    count++;
                }
            }
            if (pass == 0) found[t + 1] = count;
        }
        if (pass == 0) {
            for (int t = 0; t < threads; t++) found[t + 1] += found[t];
            if (found[threads] == 0) break;
            list_out = (radix_run_t *)ctools_safe_malloc2(found[threads], sizeof(radix_run_t));
            if (!list_out) {
                rc = STATA_ERR_MEMORY;
                break;
            }
        }
    }
    if (rc == STATA_OK && list_out) {
        *runs = list_out;
        *nruns = found[threads];
    }
    free(found);
    return rc;
}

/*
    Parallel radix_refine_run for large runs: every worker of ctx shares each
    word's gather and passes. Runs of at least `large` records that tie on a
    word are refined the same way, one after another; smaller ones are shared
    among the workers whole, one thread per run.
*/
static stata_retcode radix_refine_parallel(const radix_key_list_t *list, size_t k, size_t off,
                                           radix_records_t cur, radix_records_t tmp,
                                           perm_idx_t *out, size_t n, size_t large,
                                           radix_sort_context_t *ctx)
{
    const int threads = ctx->num_threads;

    for (;;) {
        const stata_variable *var = &list->vars[k];
        const uint64_t first = radix_first_word(cur.rows, var, off);
        uint64_t differ = 0;
        #pragma omp parallel num_threads(threads) if(threads > 1) reduction(|:differ)
        {
            int team = 1, tid = 0;
            #ifdef _OPENMP
            team = omp_get_num_threads();
            tid = omp_get_thread_num();
            #endif
            size_t start, end;
            radix_record_share(n, team, tid, &start, &end);
            differ = radix_gather_words(cur.keys + start, cur.rows + start, end - start, var, off,
                                        first);
        }
        if (radix_sort_word(&cur, &tmp, NULL, n, differ, large, ctx) == &tmp) {
            const radix_records_t swap = cur;
            cur = tmp;
            tmp = swap;
        }
        {
            const perm_idx_t *rows = cur.rows;
            #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
            for (size_t i = 0; i < n; i++) out[i] = rows[i];
        }
        if (radix_word_is_last(list, k)) return STATA_OK;

        /* One run: continue with the next word. */
        if (cur.keys[0] == cur.keys[n - 1]) {
            if (!radix_next_word(list, &k, &off, cur.keys[0])) return STATA_OK;
            continue;
        }

        radix_run_t *runs;
        size_t nruns;
        stata_retcode rc = radix_collect_runs(list, k, off, cur.keys, n, threads, &runs, &nruns);
        for (size_t r = 0; r < nruns && rc == STATA_OK; r++) {
            const size_t start = runs[r].start, len = runs[r].len;
            if (len < large) continue;
            size_t k2 = k, off2 = off;
            radix_next_word(list, &k2, &off2, cur.keys[start]);
            const radix_records_t sub = {cur.keys + start, cur.rows + start};
            const radix_records_t sub_tmp = {tmp.keys + start, tmp.rows + start};
            rc = radix_refine_parallel(list, k2, off2, sub, sub_tmp, out + start, len, large, ctx);
        }
        if (rc == STATA_OK && nruns > 0) {
            #pragma omp parallel num_threads(threads) if(threads > 1)
            {
                uint32_t counts[8][RADIX_SIZE];
                #pragma omp for schedule(dynamic, 1)
                for (size_t r = 0; r < nruns; r++) {
                    const size_t start = runs[r].start, len = runs[r].len;
                    if (len >= large) continue;
                    size_t k2 = k, off2 = off;
                    radix_next_word(list, &k2, &off2, cur.keys[start]);
                    const radix_records_t sub = {cur.keys + start, cur.rows + start};
                    const radix_records_t sub_tmp = {tmp.keys + start, tmp.rows + start};
                    radix_refine_run(list, k2, off2, sub, sub_tmp, out + start, len, counts);
                }
            }
        }
        free(runs);
        return rc;
    }
}

/*
    Stable sort of data->sort_order by columns cols[0..ncols-1] (zero-based
    positions in data->vars, most significant first), numeric or string.
*/
static stata_retcode sort_key_columns(stata_data *data, const int *cols, size_t ncols)
{
    const size_t nobs = data->nobs;
    const int threads = radix_record_threads(nobs);
    stata_retcode rc = STATA_ERR_MEMORY;
    radix_records_t a = {NULL, NULL}, b = {NULL, NULL};
    radix_sort_context_t *ctx = NULL;
    stata_variable *vars = (stata_variable *)ctools_safe_malloc2(ncols, sizeof(stata_variable));
    radix_key_list_t list = {vars, ncols};

    if (nobs < 2) {
        free(vars);
        return STATA_OK;
    }
    a.keys = (uint64_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(uint64_t));
    b.keys = (uint64_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(uint64_t));
    a.rows = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(perm_idx_t));
    b.rows = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(perm_idx_t));
    if (threads > 1) ctx = radix_context_alloc(threads);
    if (!vars || !a.keys || !b.keys || !a.rows || !b.rows || (threads > 1 && !ctx)) goto done;

    for (size_t c = 0; c < ncols; c++) {
        vars[c] = data->vars[cols[c]];
        if (vars[c].type == STATA_TYPE_DOUBLE ? !vars[c].data.dbl : !vars[c].data.str) {
            rc = STATA_ERR_INVALID_INPUT;
            goto done;
        }
    }
    memcpy(a.rows, data->sort_order, nobs * sizeof(perm_idx_t));
    if (threads > 1) {
        /* Runs and buckets this large share every worker; smaller ones take
           one each, and there are enough of them to keep the workers busy. */
        size_t large = nobs / (4 * (size_t)threads);
        if (large < (size_t)MIN_OBS_PER_THREAD * 2) large = (size_t)MIN_OBS_PER_THREAD * 2;
        rc = radix_refine_parallel(&list, 0, 0, a, b, data->sort_order, nobs, large, ctx);
    } else {
        uint32_t counts[8][RADIX_SIZE];
        radix_refine_run(&list, 0, 0, a, b, data->sort_order, nobs, counts);
        rc = STATA_OK;
    }

done:
    radix_context_free(ctx);
    ctools_aligned_free(a.keys);
    ctools_aligned_free(b.keys);
    ctools_aligned_free(a.rows);
    ctools_aligned_free(b.rows);
    free(vars);
    return rc;
}

/* Stable sort of data->sort_order by one numeric column (0-based var_idx). */
static stata_retcode sort_by_numeric_var(stata_data *data, int var_idx)
{
    return sort_key_columns(data, &var_idx, 1);
}

/*
    Sort by a single string variable using radix sort.
    Optimization 2: Pre-caches string lengths.
    Optimization 5: Uses parallel sort for large datasets.
*/
static stata_retcode sort_by_string_var(stata_data *data, int var_idx)
{
    perm_idx_t *order_a;
    perm_idx_t *order_b;
    perm_idx_t *current_order;
    perm_idx_t *temp_order;
    size_t *str_lengths;
    size_t max_len = 0;
    size_t i, len;
    int char_pos;
    char **str_data;
    int use_parallel;
    int num_threads;
    string_sort_context_t *ctx = NULL;
    int swapped = 0;

    str_data = data->vars[var_idx].data.str;

    /* Safety check for NULL string data array */
    if (str_data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Optimization 2: Pre-cache all string lengths */
    str_lengths = (size_t *)malloc(data->nobs * sizeof(size_t));
    if (str_lengths == NULL) {
        return STATA_ERR_MEMORY;
    }

    for (i = 0; i < data->nobs; i++) {
        /* Handle NULL string pointers gracefully */
        if (str_data[i] == NULL) {
            len = 0;
        } else {
            len = strlen(str_data[i]);
        }
        str_lengths[i] = len;
        if (len > max_len) {
            max_len = len;
        }
    }

    if (max_len == 0) {
        free(str_lengths);
        return STATA_OK;
    }

    /* Allocate temporary order array - aligned for cache efficiency (uses perm_idx_t) */
    order_a = data->sort_order;
    order_b = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE, data->nobs * sizeof(perm_idx_t));
    if (order_b == NULL) {
        free(str_lengths);
        return STATA_ERR_MEMORY;
    }

    /* Decide on parallelization */
    use_parallel = (data->nobs >= MIN_OBS_PER_THREAD * 2);
    num_threads = ctools_get_openmp_threads();
    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) {
            use_parallel = 0;
        }
    }

    if (use_parallel) {
        ctx = string_context_alloc(num_threads);
        if (ctx == NULL) {
            ctools_aligned_free(order_b);
            free(str_lengths);
            return STATA_ERR_MEMORY;
        }
    }

    current_order = order_a;
    temp_order = order_b;

    /* Perform radix sort passes (LSD - from rightmost character to leftmost) */
    for (char_pos = (int)max_len - 1; char_pos >= 0; char_pos--) {
        int skipped;

        if (use_parallel) {
            skipped = radix_sort_pass_string_parallel(current_order, temp_order, str_data,
                                                       str_lengths, data->nobs, (size_t)char_pos, ctx);
        } else {
            skipped = radix_sort_pass_string(current_order, temp_order, str_data,
                                             str_lengths, data->nobs, (size_t)char_pos);
        }

        if (!skipped) {
            perm_idx_t *tmp = current_order;
            current_order = temp_order;
            temp_order = tmp;
            swapped = !swapped;
        }
    }

    /* Ensure final result is in data->sort_order */
    if (current_order != data->sort_order) {
        memcpy(data->sort_order, current_order, data->nobs * sizeof(perm_idx_t));
    }

    string_context_free(ctx);
    ctools_aligned_free(order_b);
    free(str_lengths);
    return STATA_OK;
}

/*
    ctools_sort_radix_lsd_order_only - LSD radix sort without applying permutation

    Computes sort_order but does NOT apply the permutation to data.
    After this call, data->sort_order contains the permutation but data is unchanged.
    Call ctools_apply_permutation() separately to apply the permutation.

    @param data       Dataset to sort
    @param sort_vars  Array of 1-based variable indices to sort by
    @param nsort      Number of sort variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_radix_lsd_order_only(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Numeric keys and strings of up to CTOOLS_SORT_RECORD_STRING_MAX bytes
       share one record sort that stops once every tie is broken. */
    {
        int records = 1;
        int *cols = (int *)malloc(nsort * sizeof(int));
        if (cols == NULL) return STATA_ERR_MEMORY;
        for (k = 0; k < (int)nsort; k++) {
            var_idx = sort_vars[k] - 1;
            if (var_idx < 0 || var_idx >= (int)data->nvars) {
                free(cols);
                return STATA_ERR_INVALID_INPUT;
            }
            cols[k] = var_idx;
            if (data->vars[var_idx].type == STATA_TYPE_STRING &&
                data->vars[var_idx].str_maxlen > CTOOLS_SORT_RECORD_STRING_MAX) {
                records = 0;
            }
        }
        rc = records ? sort_key_columns(data, cols, nsort) : STATA_OK;
        free(cols);
        if (records) return rc;
    }

    /*
        For stable LSD radix sort with multiple keys:
        Sort from the LAST (least significant) key to the FIRST (most significant).
    */
    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;  /* Convert to 0-based index */

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = sort_by_numeric_var(data, var_idx);
        } else {
            rc = sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    /* Do NOT apply permutation - caller will do it separately */
    return STATA_OK;
}

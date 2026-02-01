/*
    csort_stream.h
    Streaming permutation for memory-efficient sorting

    This module implements a memory-efficient alternative to the standard
    load-all → sort → permute-all → store-all approach. Instead:

    1. Load ONLY key variables into C memory
    2. Sort keys to produce permutation array
    3. Apply permutation to keys in C memory (they're already loaded)
    4. Write sorted keys back to Stata
    5. Stream-permute non-key variables in blocks:
       - For each output block, gather required source rows
       - Read from Stata, permute in buffer, write back to Stata
       - Memory: O(block_size × num_non_key_vars) per block

    Memory Comparison (10M obs, 3 keys, 97 non-key vars):
    - Standard approach: ~8 GB (all data) + permutation overhead
    - Streaming approach: ~240 MB (keys only) + ~50 MB (per block buffer)

    Trade-offs:
    - Pros: Dramatically reduced memory footprint
    - Cons: Non-key variables read twice (once for each source position)
            Random access to Stata data has overhead

    Best for:
    - Large datasets where memory is constrained
    - Sorting on few key variables with many non-key variables
    - Datasets too large to fit in memory with standard approach
*/

#ifndef CSORT_STREAM_H
#define CSORT_STREAM_H

#include <stddef.h>
#include <stdint.h>
#include "ctools_types.h"

/* ============================================================================
   Configuration
   ============================================================================ */

/*
    Default block size for streaming permutation.
    Larger blocks = fewer passes, better amortization of overhead
    Smaller blocks = less memory usage

    64K observations is a good balance:
    - 64K × 100 vars × 8 bytes = ~50 MB per block
    - Small enough to fit in L3 cache for many systems
*/
#define CSORT_STREAM_DEFAULT_BLOCK_SIZE (64 * 1024)

/*
    Minimum block size to ensure efficiency.
    Below this, per-block overhead dominates.
*/
#define CSORT_STREAM_MIN_BLOCK_SIZE 4096

/*
    Threshold for using streaming vs standard approach.
    If total data size exceeds this, prefer streaming.
    Default: 1 GB (streaming helps when memory is tight)
*/
#define CSORT_STREAM_THRESHOLD_BYTES (1024ULL * 1024 * 1024)

/* ============================================================================
   Types
   ============================================================================ */

/*
    Timing breakdown for streaming sort operations.
*/
typedef struct {
    double load_keys_time;      /* Time to load key variables */
    double sort_time;           /* Time for sort algorithm */
    double permute_keys_time;   /* Time to permute key variables in C */
    double store_keys_time;     /* Time to write sorted keys to Stata */
    double stream_nonkeys_time; /* Time to stream-permute non-key variables (total) */
    /* Detailed breakdown of stream_nonkeys_time: */
    double stream_build_inv_time;   /* Time to build inverse permutation */
    double stream_scatter_time;     /* Time for Phase 1: read Stata + scatter to buffer */
    double stream_writeback_time;   /* Time for Phase 2: write buffer back to Stata */
    double stream_string_time;      /* Time for string variable processing */
    double total_time;          /* Total wall-clock time */
    size_t num_blocks;          /* Number of blocks processed */
    size_t block_size;          /* Actual block size used */
    size_t n_numeric_vars;      /* Number of numeric non-key variables */
    size_t n_string_vars;       /* Number of string non-key variables */
} csort_stream_timings;

/* ============================================================================
   Streaming Sort API
   ============================================================================ */

/*
    Perform a memory-efficient streaming sort.

    This is the main entry point for streaming sort. It:
    1. Loads only key variables
    2. Sorts to produce permutation
    3. Applies permutation to keys and writes them back
    4. Streams permutation to non-key variables in blocks

    @param key_var_indices   [in] 1-based Stata indices of key variables (sort columns)
    @param nkeys             [in] Number of key variables
    @param all_var_indices   [in] 1-based Stata indices of ALL variables in dataset
    @param nvars             [in] Total number of variables
    @param algorithm         [in] Sort algorithm to use
    @param block_size        [in] Block size for streaming (0 = use default)
    @param vars_per_batch    [in] Number of variables to process at a time (1-16)
    @param timings           [out] Optional timing breakdown (can be NULL)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode csort_stream_sort(
    int *key_var_indices,
    size_t nkeys,
    int *all_var_indices,
    size_t nvars,
    sort_algorithm_t algorithm,
    size_t block_size,
    int vars_per_batch,
    csort_stream_timings *timings
);

/*
    Choose optimal block size based on dataset characteristics.

    @param nobs              [in] Number of observations
    @param nvars_nonkey      [in] Number of non-key variables
    @param available_memory  [in] Available memory in bytes (0 = auto-detect)

    @return Recommended block size
*/
size_t csort_stream_choose_block_size(
    size_t nobs,
    size_t nvars_nonkey,
    size_t available_memory
);

/*
    Determine if streaming sort should be used over standard sort.

    @param nobs   [in] Number of observations
    @param nvars  [in] Total number of variables
    @param nkeys  [in] Number of key (sort) variables

    @return 1 if streaming is recommended, 0 otherwise
*/
int csort_stream_recommended(size_t nobs, size_t nvars, size_t nkeys);

/* ============================================================================
   Internal Functions (exposed for testing)
   ============================================================================ */

/*
    Apply permutation to non-key variables using variable-at-a-time approach.

    For each non-key variable (parallelized with OpenMP):
      1. Allocate O(nobs) buffer
      2. Read all values in permuted order: buffer[i] = Stata[perm[i]]
      3. Write all values back: Stata[i] = buffer[i]
      4. Free buffer, move to next variable

    This avoids the read-after-write hazard of block-based streaming.

    @param perm              [in] Permutation array (sorted_idx -> original_idx)
    @param nobs              [in] Number of observations
    @param nonkey_var_indices[in] 1-based Stata indices of non-key variables
    @param nvars_nonkey      [in] Number of non-key variables
    @param obs1              [in] First observation in Stata (1-based)
    @param block_size        [in] Unused (kept for API compatibility)
    @param vars_per_batch    [in] Number of variables to process at a time (1-16)
    @param timings           [out] Optional timing breakdown (can be NULL)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode csort_stream_apply_permutation(
    const perm_idx_t *perm,
    size_t nobs,
    int *nonkey_var_indices,
    size_t nvars_nonkey,
    size_t obs1,
    size_t block_size,
    int vars_per_batch,
    csort_stream_timings *timings
);

#endif /* CSORT_STREAM_H */

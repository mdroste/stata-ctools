#ifndef CTOOLS_DATA_IO_H
#define CTOOLS_DATA_IO_H
#include "ctools_data.h"

/* ---------------------------------------------------------------------------
   Data I/O Operations
   --------------------------------------------------------------------------- */

// Write all variables from C memory back to Stata
stata_retcode ctools_data_store(stata_data *data, size_t obs1);

// Write in data->sort_order without copying/permuting the loaded columns.
// Leaves both the columns and sort_order unchanged.
stata_retcode ctools_data_store_sorted(stata_data *data, size_t obs1);

// Write variables from C memory to specified Stata variable indices
// var_indices[j] is the 1-based Stata variable index for data->vars[j]
stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices, size_t nvars,
                                          size_t obs1);

// Unified store with auto-dispatch: row-parallel for nvars==1, column-parallel for nvars>=2
// var_indices == NULL → sequential 1..data->nvars
stata_retcode ctools_data_store_ex(stata_data *data, int *var_indices, size_t nvars, size_t obs1);

// Write specific output variable values to Stata using a source row mapping
// For each output row i, reads from source_rows[i] (0-based, -1 = missing)
// var_idx: 1-based Stata variable index
// Handles both numeric and string variables
stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows, size_t output_nobs,
                                         size_t obs1);

/* ---------------------------------------------------------------------------
   Filtered Data Loading

   Unified data loading with if/in filtering at load time. Eliminates per-command
   filtering code and reduces memory allocation to only filtered observations.
   --------------------------------------------------------------------------- */

/*
    Filtered dataset structure.

    Contains the standard stata_data plus observation mapping for write-back.
    Memory is allocated only for N_filtered observations, not the full in-range.

    Usage:
    1. Call ctools_data_load() to load only observations passing SF_ifobs()
    2. Access data directly: fd.data.vars[k].data.dbl[i] (no indirection needed)
    3. For write-back: ctools_store_filtered(values, N, var_idx, fd.obs_map)
    4. Call ctools_filtered_data_free() when done
*/
typedef struct
{
    stata_data data;     /* Standard data structure (N = N_filtered) */
    perm_idx_t *obs_map; /* obs_map[i] = 1-based Stata obs for filtered index i */
    size_t n_range;      /* Original range size (before filtering) */
    int was_filtered;    /* 1 if any filtering occurred, 0 if identity */
} ctools_filtered_data;

/*
    Initialize filtered data structure to safe empty state.
    Sets all pointers to NULL, all counts to 0.
*/
void ctools_filtered_data_init(ctools_filtered_data *fd);

/*
    Free all memory associated with filtered data.
    Safe to call multiple times; resets structure to empty state.
*/
void ctools_filtered_data_free(ctools_filtered_data *fd);

/*
    Flags for ctools_data_load()
*/
#define CTOOLS_LOAD_CHECK_IF 0x00      /* Default: check SF_ifobs for each observation */
#define CTOOLS_LOAD_SKIP_IF 0x01       /* Skip SF_ifobs checks, assume all observations pass */
#define CTOOLS_LOAD_NO_SORT_ORDER 0x02 /* Omit identity sort_order; caller will not sort */
#define CTOOLS_LOAD_MARKER 0x04        /* Select by a sample marker, the last plugin variable */

/*
    Load variables with if/in filtering applied at load time.

    This is the PRIMARY data loading function for ctools. It handles:
    - Selective loading of specific variables
    - "Load all" mode when var_indices is NULL
    - If/in filtering at load time (or skip with CTOOLS_LOAD_SKIP_IF)

    Performs two-pass loading:
    1. Count observations passing SF_ifobs() and build obs_map
    2. Load only filtered observations (memory = O(N_filtered × K))

    Special modes:
    - If var_indices is NULL: loads ALL variables (equivalent to old ctools_data_load)
    - If CTOOLS_LOAD_SKIP_IF flag is set: skips SF_ifobs() checks entirely
    - If CTOOLS_LOAD_MARKER is set (without SKIP_IF): the last plugin variable,
      SF_nvars(), is a numeric sample marker such as the ado's touse, passed
      in place of `if touse'. Observations of the in-range where it is nonzero
      are selected, as `if marker' selects them (missing counts as nonzero).
      The marker is read in parallel, where SF_ifobs() costs one serial SPI
      call per observation. Load-all mode leaves the marker out.
    - If CTOOLS_LOAD_NO_SORT_ORDER is set: leaves data.sort_order NULL.
      Use only when no subsequent operation requires a sort/permutation order.
      Without this flag the identity sort_order is always populated.

    @param result       [out] Filtered data structure to populate (caller frees)
    @param var_indices  [in]  Array of 1-based Stata variable indices to load,
                              or NULL to load all variables
    @param nvars        [in]  Number of variables (ignored if var_indices is NULL)
    @param obs_start    [in]  First observation (1-based), 0 = use SF_in1()
    @param obs_end      [in]  Last observation (1-based), 0 = use SF_in2()
    @param flags        [in]  CTOOLS_LOAD_CHECK_IF (0), CTOOLS_LOAD_SKIP_IF or
                              CTOOLS_LOAD_MARKER, with CTOOLS_LOAD_NO_SORT_ORDER

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if result is NULL
            STATA_ERR_MEMORY on allocation failure
*/
stata_retcode ctools_data_load(ctools_filtered_data *result, int *var_indices, size_t nvars,
                               size_t obs_start, size_t obs_end, int flags);

/*
    Resolve the if/in selection without loading variables.

    Fills result->obs_map, n_range and was_filtered, and sets
    result->data.nobs to the number of selected observations (no variables
    are allocated). Free with ctools_filtered_data_free().
*/
stata_retcode ctools_data_select(ctools_filtered_data *result, size_t obs_start, size_t obs_end,
                                 int flags);

/*
    Extended data load with string width hints for flat buffer optimization.

    Same as ctools_data_load but accepts an optional array of string variable
    widths. When provided, string variables are loaded into contiguous flat
    buffers (one SF_sdata read per observation, no arena overhead) instead
    of the default arena path (strlen + memcpy + atomic CAS per string).

    @param str_widths   [in]  Array indexed by plugin-visible variable index - 1
                              (e.g. 17 for str17, 0 for numeric/unknown/strL).
                              Array length must be SF_nvars(), including when
                              var_indices selects or reorders a subset.
                              NULL = auto-detect the caller width macro; if
                              unavailable, use the default arena path.

    All other parameters are identical to ctools_data_load.
*/
stata_retcode ctools_data_load_ex(ctools_filtered_data *result, int *var_indices, size_t nvars,
                                  size_t obs_start, size_t obs_end, int flags,
                                  const int *str_widths);

/*
    Write filtered variable values back to Stata using obs_map.

    For each filtered index i, writes values[i] to Stata observation obs_map[i].
    This enables write-back after processing without maintaining full-range arrays.

    @param values       [in] Array of values to store (length n_filtered)
    @param n_filtered   [in] Number of filtered observations
    @param var_idx      [in] 1-based Stata variable index to store to
    @param obs_map      [in] Observation mapping from filtered to Stata indices

    @return STATA_OK on success, STATA_ERR_INVALID_INPUT if parameters invalid
*/
stata_retcode ctools_store_filtered(double *values, size_t n_filtered, int var_idx,
                                    perm_idx_t *obs_map);

/*
    Row-parallel store for a single numeric variable.
    Counterpart to ctools_data_load_single_var_rowpar().
    Skips per-element bounds checks (obs_map is trusted).

    @param values       [in] Array of values to store (length n_filtered)
    @param n_filtered   [in] Number of filtered observations
    @param var_idx      [in] 1-based Stata variable index to store to
    @param obs_map      [in] Observation mapping from filtered to Stata indices
*/
stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered, int var_idx,
                                           perm_idx_t *obs_map);

/*
    Write filtered string values back to Stata using obs_map.

    @param strings      [in] Array of strings to store (length n_filtered)
    @param n_filtered   [in] Number of filtered observations
    @param var_idx      [in] 1-based Stata variable index to store to
    @param obs_map      [in] Observation mapping from filtered to Stata indices

    @return STATA_OK on success, STATA_ERR_INVALID_INPUT if parameters invalid
*/
stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered, int var_idx,
                                        perm_idx_t *obs_map);

#endif /* CTOOLS_DATA_IO_H */

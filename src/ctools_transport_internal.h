#ifndef CTOOLS_TRANSPORT_INTERNAL_H
#define CTOOLS_TRANSPORT_INTERNAL_H
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdatomic.h>
#include "stplugin.h"
#include "ctools_data_io.h"
#include "ctools_threads.h"
#include "ctools_config.h"
#include "ctools_arena.h"
/* Validated selection plan. indices are 1-based plugin columns; obs_map holds
 * 1-based Stata rows. owned_* are temporary plan storage. The receiving filtered
 * dataset takes ownership of obs_map before allocation/transfer begins. Widths
 * remain hints, never unchecked bounds. Scheduler gates stay with hot kernels. */
typedef struct
{
    int *indices, *owned_indices;
    const int *widths;
    int *owned_widths;
    size_t nvars;
    perm_idx_t *obs_map;
    size_t n_range, nobs;
    int was_filtered;
} ctools_load_plan;
/* Numeric reads use the unchecked SPI accessor: every public entry point
 * validates variable indices and observation ranges before its loops, so the
 * per-cell bounds checks were redundant (measured 1.10-1.23x on loads).
 * Define CTOOLS_IO_CHECKED_SPI to restore the checked read callback.
 *
 * Writes must stay on the checked callback: concurrent unchecked stores to
 * different variables corrupted whole columns in real Stata (cwinsor and
 * cinterpolate validation), redirecting one variable's column
 * into another, even though single-threaded semantics are bitwise identical.
 * Do not switch IO_VSTORE_FN to (_stata_)->store. */
#ifdef CTOOLS_IO_CHECKED_SPI
#define IO_VDATA_FN ((_stata_)->safevdata)
#else
#define IO_VDATA_FN ((_stata_)->vdata)
#endif
#define IO_VSTORE_FN ((_stata_)->safestore)

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045
enum
{
    IO_NUMERIC = 0,
    IO_FIXED_STRING = 1,
    IO_LONG_STRING = 2
};

/* Row tiles benefit wide fixed-string transfers, but can slow narrow/numeric
 * data. Keep the column scheduler below this conservative width/size gate. */
#define IO_MIN_TILED_VARS 4
#define IO_MIN_TILED_STRING_BYTES 256
#define IO_LOAD_TILE_ROWS 512
#define IO_STORE_TILE_ROWS 4096

/* Numeric tiling needs both enough columns and enough selected rows. Compact
 * byte columns and short, wide inputs favor column scheduling even with many
 * cells. Re-measure the tile sweep before changing these thresholds. */
#define IO_NUMERIC_TILE_MIN_VARS 128
#define IO_NUMERIC_TILE_MIN_ROWS 50000
#define IO_SERIAL_BASE_WORK 16384
#define IO_SERIAL_COLUMN_WORK 256
#define IO_DIRECT_GATHER_MAX_ROWS 262144

#define IO_SLOT_GUARD ((size_t)STATA_STR_MAXLEN + 1)
#define IO_DIRECT_READ_MIN_WIDTH 32
#ifndef IO_CHUNK_MIN_NUMERIC_ROWS
#define IO_CHUNK_MIN_NUMERIC_ROWS 1500000
#endif

typedef enum
{
    IO_MODE_LOAD,
    IO_MODE_STORE
} io_mode_t;
typedef void *(*io_thread_func)(void *);
typedef struct
{
    stata_variable *var;
    int var_idx;
    size_t obs1;           /* first Stata obs of this chunk (identity) */
    const perm_idx_t *map; /* filtered loads: 1-based source rows */
    size_t dest;           /* first destination row */
    size_t count;
    atomic_int *error;
} ctools_io_chunk;
typedef struct
{
    stata_variable *var;
    int var_idx;
    size_t obs1; /* first destination obs of this chunk */
    size_t dest; /* offset into column data / order */
    size_t count;
    size_t total_nobs;
    const perm_idx_t *order;
    atomic_int *error;
} ctools_store_chunk;
typedef struct
{
    stata_variable *var;
    int var_idx;
    perm_idx_t *obs_map;
    size_t n_filtered;
    int is_string;
    int str_width;
    int error;
} filtered_var_io_args;
stata_retcode ctools_io_validate_variable(int var_idx, int expected_string);
stata_retcode ctools_io_validate_range(size_t obs1, size_t nobs);
stata_retcode ctools_io_validate_obs_map(const perm_idx_t *obs_map, size_t nobs);
int ctools_io_read_stata_strl(ST_int var_idx, ST_int obs, char *value);
int ctools_io_variable_kind(int var_idx);
int ctools_io_read_string_slot(ST_IIIS read_fn, int var_idx, ST_int obs, char *slot, size_t width);
int ctools_io_read_string_slot_direct(ST_IIIS read_fn, int var_idx, ST_int obs, char *slot,
                                      size_t width, const char *safe_end);
stata_retcode ctools_io_load_status(int rc);
void ctools_io_record_io_error(atomic_int *error, int rc);
stata_retcode ctools_io_store_status(int rc);
int *ctools_io_ctools_read_strw_from_stata(size_t nvars);
char *ctools_io_copy_private_string(ctools_string_arena *arena, const char *value);
ctools_string_arena *ctools_io_create_private_string_arena(size_t nobs, int str_width);
int ctools_io_load_single_variable(stata_variable *var, int var_idx, size_t obs1, size_t nobs,
                                   int is_string, int str_width);
void *ctools_io_load_variable_thread(void *arg);
stata_retcode ctools_io_init_data_structure(stata_data *data, size_t nvars, size_t nobs, int flags);
void ctools_io_init_io_thread_args(ctools_var_io_args *args, stata_data *data, int *var_indices,
                                   size_t nvars, size_t obs1, size_t nobs, io_mode_t mode,
                                   const int *str_widths);
size_t ctools_io_column_work(int is_string, size_t width);
size_t ctools_io_add_work(size_t total, size_t column);
int ctools_io_serial_work(size_t nobs, size_t nvars, size_t row_work);
int ctools_io_execute_io_parallel(ctools_var_io_args *args, size_t nvars, io_thread_func func);
void *ctools_io_load_chunk_thread(void *arg);
size_t ctools_io_chunks_per_column(size_t nvars, size_t nobs, size_t min_chunk_rows);
int ctools_io_try_chunked_load(stata_data *data, int *var_indices, size_t nvars,
                               const perm_idx_t *obs_map, size_t nobs, int was_filtered,
                               stata_retcode *rc_out);
int ctools_io_store_single_variable(stata_variable *var, int var_idx, size_t obs1, size_t nobs,
                                    const perm_idx_t *order);
void *ctools_io_store_variable_thread(void *arg);
void *ctools_io_store_chunk_thread(void *arg);
stata_retcode ctools_io_data_store_impl(stata_data *data, int *var_indices, size_t nvars,
                                        size_t obs1, const perm_idx_t *order);
stata_retcode ctools_data_store_ex(stata_data *data, int *var_indices, size_t nvars, size_t obs1);
stata_retcode ctools_data_store_sorted(stata_data *data, size_t obs1);
stata_retcode ctools_data_store(stata_data *data, size_t obs1);
stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices, size_t nvars,
                                          size_t obs1);
stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows, size_t output_nobs,
                                         size_t obs1);
stata_retcode ctools_io_build_obs_map(size_t obs_start, size_t obs_end, int flags,
                                      perm_idx_t **obs_map_out, size_t *n_filtered_out,
                                      size_t *n_range_out, int *was_filtered_out);
int ctools_io_load_row_parallel(stata_variable *var, int var_idx, perm_idx_t *obs_map,
                                size_t n_filtered, int str_width);
void ctools_filtered_data_init(ctools_filtered_data *fd);
void ctools_filtered_data_free(ctools_filtered_data *fd);
int ctools_io_probe_for_identity(ST_int obs1, ST_int obs2, size_t n_range);
int ctools_io_load_filtered_variable(stata_variable *var, int var_idx, perm_idx_t *obs_map,
                                     size_t n_filtered, int is_string, int str_width);
void *ctools_io_load_filtered_variable_thread(void *arg);
int ctools_io_load_tiled_variables(stata_data *data, int *indices, const perm_idx_t *map,
                                   const int *widths);
stata_retcode ctools_data_select(ctools_filtered_data *result, size_t obs_start, size_t obs_end,
                                 int flags);
stata_retcode ctools_data_load_ex(ctools_filtered_data *result, int *var_indices, size_t nvars,
                                  size_t obs_start, size_t obs_end, int flags,
                                  const int *str_widths);
stata_retcode ctools_data_load(ctools_filtered_data *result, int *var_indices, size_t nvars,
                               size_t obs_start, size_t obs_end, int flags);
stata_retcode ctools_store_filtered(double *values, size_t n_filtered, int var_idx,
                                    perm_idx_t *obs_map);
stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered, int var_idx,
                                           perm_idx_t *obs_map);
stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered, int var_idx,
                                        perm_idx_t *obs_map);
stata_retcode ctools_io_plan_load(ctools_load_plan *plan, int *var_indices, size_t nvars,
                                  size_t obs_start, size_t obs_end, int flags,
                                  const int *str_widths);
#endif

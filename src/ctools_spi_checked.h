/*
 * ctools_spi_checked.h - Error-checking wrappers for Stata Plugin Interface
 *
 * All SF_* write functions return ST_int error codes that are typically
 * ignored. These wrappers check return values and report errors, preventing
 * silent failures that can corrupt Stata's memory.
 *
 * Usage:
 *   #include "ctools_spi_checked.h"
 *
 *   // Instead of: SF_scal_save("__myvar", 123.0);
 *   // Use:        ctools_scal_save("__myvar", 123.0);
 *
 * Return values:
 *   0 = success
 *   non-zero = Stata error code (e.g., matrix not found, invalid index)
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_SPI_CHECKED_H
#define CTOOLS_SPI_CHECKED_H

#include "stplugin.h"

/*
 * Save a scalar value to Stata with error checking.
 * Returns 0 on success, non-zero on error.
 */
ST_retcode ctools_scal_save(const char *name, ST_double val);

/*
 * Store a value in a Stata matrix with error checking.
 * Row and column indices are 1-based (Stata convention).
 * Returns 0 on success, non-zero on error.
 */
ST_retcode ctools_mat_store(const char *name, ST_int row, ST_int col, ST_double val);

/*
 * Save a macro value to Stata with error checking.
 * Returns 0 on success, non-zero on error.
 */
ST_retcode ctools_macro_save(const char *name, const char *val);

/*
 * NOTE: No wrappers for SF_vstore/SF_sstore - these are called O(N) times
 * in tight loops where per-call overhead is unacceptable. Use SF_vstore
 * and SF_sstore directly for variable data storage.
 */

/*
 * Batch store values to a Stata matrix with error checking.
 * Stores values from a contiguous array into matrix rows/columns.
 * Returns 0 on success, first error code encountered on failure.
 *
 * Parameters:
 *   name   - Matrix name in Stata
 *   values - Contiguous array of values (row-major order)
 *   nrows  - Number of rows
 *   ncols  - Number of columns
 */
ST_retcode ctools_mat_store_batch(const char *name, const ST_double *values,
                                   ST_int nrows, ST_int ncols);

/*
 * Store a vector (1D array) as a row vector in a Stata matrix.
 * Returns 0 on success, first error code encountered on failure.
 */
ST_retcode ctools_mat_store_rowvec(const char *name, const ST_double *values, ST_int n);

/*
 * Store a vector (1D array) as a column vector in a Stata matrix.
 * Returns 0 on success, first error code encountered on failure.
 */
ST_retcode ctools_mat_store_colvec(const char *name, const ST_double *values, ST_int n);

#endif /* CTOOLS_SPI_CHECKED_H */

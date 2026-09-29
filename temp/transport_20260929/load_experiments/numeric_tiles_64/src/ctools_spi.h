/*
 * ctools_spi.h - Stata Plugin Interface (SPI) helpers for ctools
 *
 * Provides:
 * - Error-checking wrappers for scalar/matrix/macro storage
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_SPI_H
#define CTOOLS_SPI_H

#include "stplugin.h"

/* ============================================================================
   Error-Checking Wrappers
   ============================================================================ */

/*
 * Save a scalar value to Stata with error checking.
 * Returns 0 on success, non-zero on error.
 */
ST_retcode ctools_scal_save(const char *name, ST_double val);

/*
 * Store a value in a Stata matrix with bounds checking, including SD_FASTMODE.
 * Row and column indices are 1-based (Stata convention).
 * Returns 0 on success, non-zero on error.
 */
ST_retcode ctools_mat_store(const char *name, ST_int row, ST_int col, ST_double val);

/*
 * NOTE: No wrappers for SF_vstore/SF_sstore - these are called O(N) times
 * in tight loops where per-call overhead is unacceptable. Use SF_vstore
 * and SF_sstore directly for variable data storage.
 */

#endif /* CTOOLS_SPI_H */

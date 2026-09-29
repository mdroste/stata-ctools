/*
 * ctools_spi.c - Error-checking wrappers for Stata Plugin Interface
 *
 * Provides wrappers around SF_* functions that check return codes and
 * report errors. This helps prevent silent failures that can lead to
 * memory corruption in Stata.
 *
 * Part of the ctools suite for Stata.
 */

#include "ctools_spi.h"
#include <stdio.h>

/*
 * Save a scalar value to Stata with error checking.
 */
ST_retcode ctools_scal_save(const char *name, ST_double val)
{
    ST_retcode rc = SF_scal_save((char *)name, val);
    if (rc != 0) {
        char buf[128];
        snprintf(buf, sizeof(buf), "ctools: Failed to save scalar '%s' (rc=%d)\n",
                 name ? name : "(null)", (int)rc);
        SF_error(buf);
    }
    return rc;
}

/*
 * Store a value in a Stata matrix with error checking.
 */
ST_retcode ctools_mat_store(const char *name, ST_int row, ST_int col, ST_double val)
{
    /* SF_mat_store selects the unchecked callback under SD_FASTMODE.
     * Result matrices belong to Stata: a bad shape must return an error,
     * never write beyond the host allocation. */
    ST_retcode rc = (_stata_)->safematstore((char *)name, row, col, val);
    if (rc != 0) {
        char buf[128];
        snprintf(buf, sizeof(buf), "ctools: Failed to store matrix '%s'[%d,%d] (rc=%d)\n",
                 name ? name : "(null)", (int)row, (int)col, (int)rc);
        SF_error(buf);
    }
    return rc;
}

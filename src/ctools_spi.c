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
    ST_retcode rc = SF_mat_store((char *)name, row, col, val);
    if (rc != 0) {
        char buf[128];
        snprintf(buf, sizeof(buf), "ctools: Failed to store matrix '%s'[%d,%d] (rc=%d)\n",
                 name ? name : "(null)", (int)row, (int)col, (int)rc);
        SF_error(buf);
    }
    return rc;
}

/*
 * Save a macro value to Stata with error checking.
 */
ST_retcode ctools_macro_save(const char *name, const char *val)
{
    ST_retcode rc = SF_macro_save((char *)name, (char *)val);
    if (rc != 0) {
        char buf[128];
        snprintf(buf, sizeof(buf), "ctools: Failed to save macro '%s' (rc=%d)\n",
                 name ? name : "(null)", (int)rc);
        SF_error(buf);
    }
    return rc;
}

/*
 * Batch store values to a Stata matrix with error checking.
 * Values are stored in row-major order (C convention).
 */
ST_retcode ctools_mat_store_batch(const char *name, const ST_double *values,
                                   ST_int nrows, ST_int ncols)
{
    if (name == NULL || values == NULL) {
        SF_error("ctools: NULL pointer passed to ctools_mat_store_batch\n");
        return 198;
    }

    ST_retcode rc;
    for (ST_int i = 0; i < nrows; i++) {
        for (ST_int j = 0; j < ncols; j++) {
            rc = SF_mat_store((char *)name, i + 1, j + 1, values[i * ncols + j]);
            if (rc != 0) {
                char buf[128];
                snprintf(buf, sizeof(buf),
                         "ctools: Failed to store matrix '%s'[%d,%d] (rc=%d)\n",
                         name, (int)(i + 1), (int)(j + 1), (int)rc);
                SF_error(buf);
                return rc;
            }
        }
    }
    return 0;
}

/*
 * Store a vector as a row vector (1 x n matrix) in Stata.
 */
ST_retcode ctools_mat_store_rowvec(const char *name, const ST_double *values, ST_int n)
{
    if (name == NULL || values == NULL) {
        SF_error("ctools: NULL pointer passed to ctools_mat_store_rowvec\n");
        return 198;
    }

    ST_retcode rc;
    for (ST_int j = 0; j < n; j++) {
        rc = SF_mat_store((char *)name, 1, j + 1, values[j]);
        if (rc != 0) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: Failed to store row vector '%s'[1,%d] (rc=%d)\n",
                     name, (int)(j + 1), (int)rc);
            SF_error(buf);
            return rc;
        }
    }
    return 0;
}

/*
 * Store a vector as a column vector (n x 1 matrix) in Stata.
 */
ST_retcode ctools_mat_store_colvec(const char *name, const ST_double *values, ST_int n)
{
    if (name == NULL || values == NULL) {
        SF_error("ctools: NULL pointer passed to ctools_mat_store_colvec\n");
        return 198;
    }

    ST_retcode rc;
    for (ST_int i = 0; i < n; i++) {
        rc = SF_mat_store((char *)name, i + 1, 1, values[i]);
        if (rc != 0) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: Failed to store column vector '%s'[%d,1] (rc=%d)\n",
                     name, (int)(i + 1), (int)rc);
            SF_error(buf);
            return rc;
        }
    }
    return 0;
}

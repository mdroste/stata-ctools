/*
 * ctools_spi.h - Stata Plugin Interface (SPI) batch access macros
 *
 * Provides optimized batch read/write macros for Stata data access.
 * These macros reduce function call overhead through loop unrolling.
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_SPI_H
#define CTOOLS_SPI_H

#include "stplugin.h"

/* ============================================================================
   Batched Numeric Read Macros
   ============================================================================ */

/*
 * Batch read 8 numeric values from Stata.
 * Reduces function call overhead through macro expansion.
 */
#define SF_VDATA_BATCH8(var, obs_base, out) do { \
    SF_vdata((var), (ST_int)(obs_base),       &(out)[0]); \
    SF_vdata((var), (ST_int)((obs_base) + 1), &(out)[1]); \
    SF_vdata((var), (ST_int)((obs_base) + 2), &(out)[2]); \
    SF_vdata((var), (ST_int)((obs_base) + 3), &(out)[3]); \
    SF_vdata((var), (ST_int)((obs_base) + 4), &(out)[4]); \
    SF_vdata((var), (ST_int)((obs_base) + 5), &(out)[5]); \
    SF_vdata((var), (ST_int)((obs_base) + 6), &(out)[6]); \
    SF_vdata((var), (ST_int)((obs_base) + 7), &(out)[7]); \
} while(0)

/*
 * Batch read 16 numeric values from Stata (aggressive unrolling).
 */
#define SF_VDATA_BATCH16(var, obs_base, out) do { \
    SF_vdata((var), (ST_int)(obs_base),        &(out)[0]); \
    SF_vdata((var), (ST_int)((obs_base) + 1),  &(out)[1]); \
    SF_vdata((var), (ST_int)((obs_base) + 2),  &(out)[2]); \
    SF_vdata((var), (ST_int)((obs_base) + 3),  &(out)[3]); \
    SF_vdata((var), (ST_int)((obs_base) + 4),  &(out)[4]); \
    SF_vdata((var), (ST_int)((obs_base) + 5),  &(out)[5]); \
    SF_vdata((var), (ST_int)((obs_base) + 6),  &(out)[6]); \
    SF_vdata((var), (ST_int)((obs_base) + 7),  &(out)[7]); \
    SF_vdata((var), (ST_int)((obs_base) + 8),  &(out)[8]); \
    SF_vdata((var), (ST_int)((obs_base) + 9),  &(out)[9]); \
    SF_vdata((var), (ST_int)((obs_base) + 10), &(out)[10]); \
    SF_vdata((var), (ST_int)((obs_base) + 11), &(out)[11]); \
    SF_vdata((var), (ST_int)((obs_base) + 12), &(out)[12]); \
    SF_vdata((var), (ST_int)((obs_base) + 13), &(out)[13]); \
    SF_vdata((var), (ST_int)((obs_base) + 14), &(out)[14]); \
    SF_vdata((var), (ST_int)((obs_base) + 15), &(out)[15]); \
} while(0)

/* ============================================================================
   Batched Numeric Write Macros
   ============================================================================ */

/*
 * Batch write 8 numeric values to Stata.
 */
#define SF_VSTORE_BATCH8(var, obs_base, in) do { \
    SF_vstore((var), (ST_int)(obs_base),       (in)[0]); \
    SF_vstore((var), (ST_int)((obs_base) + 1), (in)[1]); \
    SF_vstore((var), (ST_int)((obs_base) + 2), (in)[2]); \
    SF_vstore((var), (ST_int)((obs_base) + 3), (in)[3]); \
    SF_vstore((var), (ST_int)((obs_base) + 4), (in)[4]); \
    SF_vstore((var), (ST_int)((obs_base) + 5), (in)[5]); \
    SF_vstore((var), (ST_int)((obs_base) + 6), (in)[6]); \
    SF_vstore((var), (ST_int)((obs_base) + 7), (in)[7]); \
} while(0)

/*
 * Batch write 16 numeric values to Stata (aggressive unrolling).
 */
#define SF_VSTORE_BATCH16(var, obs_base, in) do { \
    SF_vstore((var), (ST_int)(obs_base),        (in)[0]); \
    SF_vstore((var), (ST_int)((obs_base) + 1),  (in)[1]); \
    SF_vstore((var), (ST_int)((obs_base) + 2),  (in)[2]); \
    SF_vstore((var), (ST_int)((obs_base) + 3),  (in)[3]); \
    SF_vstore((var), (ST_int)((obs_base) + 4),  (in)[4]); \
    SF_vstore((var), (ST_int)((obs_base) + 5),  (in)[5]); \
    SF_vstore((var), (ST_int)((obs_base) + 6),  (in)[6]); \
    SF_vstore((var), (ST_int)((obs_base) + 7),  (in)[7]); \
    SF_vstore((var), (ST_int)((obs_base) + 8),  (in)[8]); \
    SF_vstore((var), (ST_int)((obs_base) + 9),  (in)[9]); \
    SF_vstore((var), (ST_int)((obs_base) + 10), (in)[10]); \
    SF_vstore((var), (ST_int)((obs_base) + 11), (in)[11]); \
    SF_vstore((var), (ST_int)((obs_base) + 12), (in)[12]); \
    SF_vstore((var), (ST_int)((obs_base) + 13), (in)[13]); \
    SF_vstore((var), (ST_int)((obs_base) + 14), (in)[14]); \
    SF_vstore((var), (ST_int)((obs_base) + 15), (in)[15]); \
} while(0)

#endif /* CTOOLS_SPI_H */

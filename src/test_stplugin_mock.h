/*
 * test_stplugin_mock.h
 * Mock Stata Plugin Interface for standalone testing with AddressSanitizer
 */

#ifndef TEST_STPLUGIN_MOCK_H
#define TEST_STPLUGIN_MOCK_H

#include <math.h>

/* Stata types */
typedef signed char    ST_sbyte;
typedef unsigned char  ST_ubyte;
typedef int            ST_int;
typedef unsigned       ST_unsigned;
typedef short int      ST_int2;
typedef int            ST_int4;
typedef long           ST_long;
typedef unsigned int   ST_uint4;
typedef float          ST_float;
typedef double         ST_double;
typedef unsigned char  ST_boolean;
typedef int            ST_retcode;
typedef double *       ST_dmkey;

#define bTrue  1
#define bFalse 0

/* STATA_OK is defined as enum in ctools_types.h, don't redefine */

/* Missing value check */
#define SF_is_missing(x) (isnan(x) || (x) > 8.988e307)
#define SV_missval (8.988e307 + 1)

/* Matrix info struct */
typedef struct {
    ST_int type;
    ST_int nel;
    ST_int m;
    ST_int n;
} ST_matinfo;

/* Mock SPI function declarations */
ST_int SF_in1(void);
ST_int SF_in2(void);
ST_int SF_nvars(void);
ST_int SF_nobs(void);
ST_int SF_ifobs(ST_int obs);

ST_int SF_scal_use(const char *name, ST_double *value);
ST_int SF_scal_save(const char *name, ST_double value);

ST_int SF_vdata(ST_int var, ST_int obs, ST_double *value);
ST_int SF_vstore(ST_int var, ST_int obs, ST_double value);

ST_int SF_sdata(ST_int var, ST_int obs, char *value);
ST_int SF_sstore(ST_int var, ST_int obs, char *value);

ST_int SF_var_is_string(ST_int var);
ST_int SF_var_is_strl(ST_int var);

ST_int SF_mat_store(const char *name, ST_int row, ST_int col, ST_double value);
ST_int SF_mat_el(const char *name, ST_int row, ST_int col, ST_double *value);

void SF_error(const char *msg);
void SF_display(const char *msg);

/* Macro for STDLL */
#define STDLL ST_retcode

#endif /* TEST_STPLUGIN_MOCK_H */

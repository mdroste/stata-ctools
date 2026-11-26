/*
    cegen_store.c
    Storage types and the store phase of cegen

    A result is computed in double precision and rounded to the new
    variable's storage type here, the way `generate type newvar = ...' rounds
    it. The store writes every observation of the dataset in one parallel
    pass (the result where the sample has one, the fill value elsewhere), so
    the ado can create the variable without initial values.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_config.h"
#include "cegen_internal.h"

/* Largest values of Stata's float and double storage types */
#define CEGEN_FLOAT_MAX 0x1.fffffep+126
#define CEGEN_DOUBLE_MAX 0x1.fffffffffffffp+1022

static const struct
{
    const char *name;
    cegen_type type;
} cegen_type_names[] = {{"byte", CEGEN_T_BYTE},   {"int", CEGEN_T_INT},       {"long", CEGEN_T_LONG},
                        {"float", CEGEN_T_FLOAT}, {"double", CEGEN_T_DOUBLE}, {"str", CEGEN_T_STR}};

int cegen_type_parse(const char *name, cegen_type *type)
{
    for (size_t i = 0; i < sizeof(cegen_type_names) / sizeof(cegen_type_names[0]); i++)
    {
        if (!strcmp(name, cegen_type_names[i].name))
        {
            *type = cegen_type_names[i].type;
            return 1;
        }
    }
    return 0;
}

const char *cegen_type_name(cegen_type type)
{
    for (size_t i = 0; i < sizeof(cegen_type_names) / sizeof(cegen_type_names[0]); i++)
        if (cegen_type_names[i].type == type)
            return cegen_type_names[i].name;
    return "";
}

static inline double cegen_round_to(double value, cegen_type type, double missval)
{
    if (!(value < missval))
        return value >= missval ? value : missval; /* missing codes; NaN */
    if (value < -CEGEN_DOUBLE_MAX)
        return missval;
    switch (type)
    {
    case CEGEN_T_BYTE:
        value = trunc(value);
        return value >= -127.0 && value <= 100.0 ? value : missval;
    case CEGEN_T_INT:
        value = trunc(value);
        return value >= -32767.0 && value <= 32740.0 ? value : missval;
    case CEGEN_T_LONG:
        value = trunc(value);
        return value >= -2147483647.0 && value <= 2147483620.0 ? value : missval;
    case CEGEN_T_FLOAT:
    {
        double rounded = (double)(float)value;
        return fabs(rounded) <= CEGEN_FLOAT_MAX ? rounded : missval;
    }
    default:
        return value;
    }
}

double cegen_type_round(double value, cegen_type type)
{
    return cegen_round_to(value, type, SV_missval);
}

void cegen_result_range(const double *values, size_t n, cegen_range *range)
{
    const double missval = SV_missval;
    size_t nonmissing = 0;
    int fractional = 0, wide = 0;
    double lo = INFINITY, hi = -INFINITY;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) \
    reduction(+ : nonmissing) reduction(| : fractional, wide) reduction(min : lo) reduction(max : hi)
    for (size_t i = 0; i < n; i++)
    {
        const double v = values[i];
        if (!(v < missval))
            continue;
        nonmissing++;
        fractional |= (v != trunc(v));
        wide |= ((double)(float)v != v);
        if (v < lo)
            lo = v;
        if (v > hi)
            hi = v;
    }
    range->nonmissing = nonmissing;
    range->all_integer = !fractional;
    range->all_float = !wide;
    range->min = nonmissing ? lo : 0.0;
    range->max = nonmissing ? hi : 0.0;
}

/* obs_map must name increasing observations of 1..nobs_total */
static int cegen_map_valid(const perm_idx_t *obs_map, size_t n, size_t nobs_total)
{
    if (!n)
        return 1;
    if (!obs_map)
        return 0;
    unsigned bad = (obs_map[0] < 1) | ((size_t)obs_map[n - 1] > nobs_total);
    for (size_t i = 1; i < n; i++)
        bad |= (obs_map[i] <= obs_map[i - 1]);
    return !bad;
}

/* The first row whose observation is at or after obs */
static size_t cegen_map_lower_bound(const perm_idx_t *obs_map, size_t n, size_t obs)
{
    size_t lo = 0, hi = n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if ((size_t)obs_map[mid] < obs)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

static stata_retcode cegen_store_status(int io_rc)
{
    if (!io_rc)
        return STATA_OK;
    return io_rc == 1 ? STATA_ERR_CANCELLED : STATA_ERR_STATA_WRITE;
}

stata_retcode cegen_store_numeric(const double *values, size_t n, const perm_idx_t *obs_map,
                                  int var, size_t nobs_total, cegen_type type, double fill,
                                  size_t *nmissing)
{
    *nmissing = 0;
    if (var < 1 || var > SF_nvars() || SF_var_is_string(var) || nobs_total > (size_t)SF_nobs() ||
        (n && !values) || !cegen_map_valid(obs_map, n, nobs_total))
        return STATA_ERR_INVALID_INPUT;
    if (!nobs_total)
        return STATA_OK;

    const double missval = SV_missval;
    const ST_IIID store_fn = (_stata_)->safestore;
    const int threads = cegen_threads_for(nobs_total);
    size_t missing = 0;
    int error = 0;
#pragma omp parallel num_threads(threads) if (threads > 1) reduction(+ : missing)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        /* Each worker owns a contiguous range of observations */
        const size_t base = nobs_total / team, extra = nobs_total % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t stop = start + base + (tid < extra);
        size_t row = cegen_map_lower_bound(obs_map, n, start + 1);
        for (size_t obs = start + 1; obs <= stop; obs++)
        {
            double v = fill;
            if (row < n && (size_t)obs_map[row] == obs)
                v = cegen_round_to(values[row++], type, missval);
            missing += !(v < missval);
            int io_rc = store_fn(var, (ST_int)obs, v);
            if (io_rc)
            {
#pragma omp atomic write
                error = io_rc;
                break;
            }
        }
    }
    *nmissing = missing;
    return cegen_store_status(error);
}

stata_retcode cegen_store_strings(char *const *values, size_t n, const perm_idx_t *obs_map,
                                  int var, size_t nobs_total, size_t *nmissing)
{
    *nmissing = 0;
    if (var < 1 || var > SF_nvars() || !SF_var_is_string(var) || nobs_total > (size_t)SF_nobs() ||
        (n && !values) || !cegen_map_valid(obs_map, n, nobs_total))
        return STATA_ERR_INVALID_INPUT;
    if (!nobs_total)
        return STATA_OK;

    const ST_IIIS store_fn = (_stata_)->sstore;
    const int threads = cegen_threads_for(nobs_total);
    size_t missing = 0;
    int error = 0;
    static char empty[1] = "";
#pragma omp parallel num_threads(threads) if (threads > 1) reduction(+ : missing)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = nobs_total / team, extra = nobs_total % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t stop = start + base + (tid < extra);
        size_t row = cegen_map_lower_bound(obs_map, n, start + 1);
        for (size_t obs = start + 1; obs <= stop; obs++)
        {
            char *text = empty;
            if (row < n && (size_t)obs_map[row] == obs)
            {
                if (values[row])
                    text = values[row];
                row++;
            }
            missing += (text[0] == '\0');
            int io_rc = store_fn(var, (ST_int)obs, text);
            if (io_rc)
            {
#pragma omp atomic write
                error = io_rc;
                break;
            }
        }
    }
    *nmissing = missing;
    return cegen_store_status(error);
}

/*
    cegen_internal.h
    Private interfaces shared by the cegen modules

    cegen_impl.c     request parsing, the compute and store phases
    cegen_groups.c   by-groups: runs of equal keys, sorted only when needed
    cegen_stats.c    statistics of one source within groups
    cegen_order.c    functions of the sorted values: rank(), mode()
    cegen_rows.c     functions across the variables of one observation
    cegen_keys.c     group(), tag(), cut(), fill(), seq(), ma()
    cegen_strings.c  concat(), ends()
    cegen_store.c    storage types and the store phase

    Indices: rows are zero-based positions in the loaded sample (load order
    is observation order); obs_map gives each row's one-based observation.
    Columns of the loaded data are zero based: sources first, then the
    by-variables.

    License: MIT
*/

#ifndef CEGEN_INTERNAL_H
#define CEGEN_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "stplugin.h"
#include "ctools_types.h"

#define CEGEN_MODULE "cegen"

/* ---------------------------------------------------------------------------
   Requests
   --------------------------------------------------------------------------- */

typedef enum
{
    CEGEN_FN_NONE = 0,
    /* statistics of one source within by-groups (cegen_stats.c) */
    CEGEN_FN_COUNT,
    CEGEN_FN_MEAN,
    CEGEN_FN_TOTAL,
    CEGEN_FN_SD,
    CEGEN_FN_MIN,
    CEGEN_FN_MAX,
    CEGEN_FN_MDEV,
    CEGEN_FN_SKEW,
    CEGEN_FN_KURT,
    CEGEN_FN_STD,
    CEGEN_FN_PC,
    CEGEN_FN_MEDIAN,
    CEGEN_FN_PCTILE,
    CEGEN_FN_IQR,
    CEGEN_FN_MAD,
    CEGEN_FN_SEQ,
    /* functions of the sorted values within by-groups (cegen_order.c) */
    CEGEN_FN_RANK,
    CEGEN_FN_MODE,
    /* functions across variables (cegen_rows.c) */
    CEGEN_FN_ROWTOTAL,
    CEGEN_FN_ROWMEAN,
    CEGEN_FN_ROWSD,
    CEGEN_FN_ROWMIN,
    CEGEN_FN_ROWMAX,
    CEGEN_FN_ROWFIRST,
    CEGEN_FN_ROWLAST,
    CEGEN_FN_ROWMISS,
    CEGEN_FN_ROWNONMISS,
    CEGEN_FN_ROWMEDIAN,
    CEGEN_FN_ROWPCTILE,
    CEGEN_FN_ANYCOUNT,
    CEGEN_FN_ANYMATCH,
    CEGEN_FN_ANYVALUE,
    CEGEN_FN_DIFF,
    /* keys and sequences (cegen_keys.c) */
    CEGEN_FN_GROUP,
    CEGEN_FN_TAG,
    CEGEN_FN_CUT,
    CEGEN_FN_FILL,
    CEGEN_FN_MA,
    /* strings (cegen_strings.c) */
    CEGEN_FN_CONCAT,
    CEGEN_FN_ENDS
} cegen_fn;

/* Stata storage types, narrowest first */
typedef enum
{
    CEGEN_T_BYTE = 0,
    CEGEN_T_INT,
    CEGEN_T_LONG,
    CEGEN_T_FLOAT,
    CEGEN_T_DOUBLE,
    CEGEN_T_STR
} cegen_type;

typedef struct
{
    cegen_fn fn;
    size_t nsrc; /* source columns: 0..nsrc-1 */
    size_t nby;  /* by columns: nsrc..nsrc+nby-1 */
    cegen_type type; /* storage type asked for (stepwise rounding uses it) */
    int xfloat;      /* working copies egen makes are float (c(type) is float) */
    int stats;       /* report the range and kind of the numeric result */

    /* function options */
    int missing;  /* total, min, max, mode, group, tag, rowtotal */
    int prop;     /* pc */
    int field, track, unique; /* rank */
    int minmode, maxmode;     /* mode */
    long long nummode;        /* mode: 0 when not given */
    int strok;    /* rownonmiss */
    int icodes;   /* cut */
    int label;    /* group: report one observation of each group */
    int nomiss;   /* ma */
    int head, last, tail, trim; /* ends */
    double p;     /* pctile, rowpctile */
    double std_mean, std_sd; /* std */
    double seq_from, seq_to, seq_block; /* seq */
    long long nq; /* cut: group(#), 0 with at() */
    long long window; /* ma: t() */
    double nobs_total; /* _N */

    double *list; /* values(), at(), fill(): owned */
    size_t nlist;
    char *punct; /* concat, ends: owned, never NULL after parsing */
    size_t punct_len;
} cegen_job;

/* What the compute phase tells the ado besides the result itself */
typedef struct
{
    size_t ngroups;
    size_t strlen_max; /* string result: bytes the new variable needs */
    /* mode() */
    int all_missing;   /* every value of the sample is missing */
    int group_missing; /* some group has only missing values */
    int multiple;      /* some group has no single mode */
    size_t max_modes;  /* most modes in a group */
    size_t bad_groups; /* groups with fewer than nummode() modes */
    /* concat(): numeric sources Stata must format, as a 0/1 text per source */
    char *needfmt;
    /* cut(): the cut points used (owned) */
    double *cuts;
    size_t ncuts;
    /* group(): observation of the first row of each group (owned) */
    perm_idx_t *reps;
    size_t nreps;
    int reps_overflow; /* more groups than a value label can hold */
    /* the storage type the result must have, when the kernel determined it */
    int has_final_type;
    cegen_type final_type;
} cegen_report;

/* The result of the compute phase: one value per loaded row */
typedef struct
{
    double *num;  /* numeric result (aligned), or NULL */
    char **str;   /* string result: one pointer per row, or NULL */
    char *arena;  /* owned text the pointers refer to, or NULL (then they
                     refer to the loaded sources, which must be kept) */
    double fill;  /* value of the observations outside the sample */
} cegen_result;

/* ---------------------------------------------------------------------------
   Groups (cegen_groups.c)
   --------------------------------------------------------------------------- */

/*
    The by-groups of the loaded rows, in one or both of two forms.

    Runs: group g holds the rows at positions starts[g] .. starts[g+1]-1 of
    order, in observation order within the group. order is NULL when the rows
    are already grouped in load order (no by-variables, or data sorted by
    them), and position is then the row.

    Numbers: ids[row] is the group of the row, 1 .. ngroups. Set when the
    groups were numbered without sorting (numeric keys of unsorted data);
    the numbering is then arbitrary, not key order.
*/
typedef struct
{
    size_t ngroups;
    size_t *starts;          /* ngroups + 1 positions (owned), or NULL */
    const perm_idx_t *order; /* rows in group order, or NULL */
    perm_idx_t *owned_order; /* order when it was built here (owned) */
    ST_int *ids;             /* group of each row (owned), or NULL */
    size_t max_count;        /* rows of the largest group (with runs) */
} cegen_groups;

/* cegen_groups_build flags */
#define CEGEN_GROUPS_RUNS 0x01      /* the caller reads starts and order */
#define CEGEN_GROUPS_KEY_ORDER 0x02 /* groups must ascend in key order */

/*
    Group the loaded rows by columns first_key .. first_key+nkeys-1 (exact
    values: distinct missing codes are distinct groups, strings compare by
    their bytes). nkeys 0 gives one group.

    Rows already in nondecreasing key order are grouped by one pass over
    neighbours. Otherwise numeric keys are numbered through the shared
    parallel remap (no sort); runs, when asked for, then come from one
    stable counting pass. String keys, and CEGEN_GROUPS_KEY_ORDER, sort
    data->sort_order instead. Without CEGEN_GROUPS_RUNS the result may hold
    only ids.
*/
stata_retcode cegen_groups_build(stata_data *data, size_t first_key, size_t nkeys, int flags,
                                 cegen_groups *groups);
void cegen_groups_free(cegen_groups *groups);

/*
    Cut the groups into at most max_chunks runs of consecutive groups with
    about the same number of rows. Returns the number of chunks; chunk c holds
    groups cuts[c] .. cuts[c+1]-1. cuts must hold max_chunks + 1 entries.
*/
size_t cegen_groups_chunks(const cegen_groups *groups, size_t max_chunks, size_t *cuts);

/* 1 when the rows are in nondecreasing order of the nkeys key columns */
int cegen_rows_sorted(const stata_data *data, const int *keys, size_t nkeys);

/* ---------------------------------------------------------------------------
   Kernels: each fills result->num or result->str for every loaded row and
   returns STATA_OK or an error. data holds the loaded columns.
   --------------------------------------------------------------------------- */

stata_retcode cegen_stats_run(const cegen_job *job, stata_data *data, cegen_result *result,
                              cegen_report *report);
stata_retcode cegen_order_run(const cegen_job *job, stata_data *data, cegen_result *result,
                              cegen_report *report);
stata_retcode cegen_rows_run(const cegen_job *job, stata_data *data, cegen_result *result,
                             cegen_report *report);
stata_retcode cegen_keys_run(const cegen_job *job, stata_data *data, const perm_idx_t *obs_map,
                             cegen_result *result, cegen_report *report);
stata_retcode cegen_strings_run(const cegen_job *job, stata_data *data, cegen_result *result,
                                cegen_report *report);

/* ---------------------------------------------------------------------------
   Shared numeric helpers (cegen_stats.c)
   --------------------------------------------------------------------------- */

/*
    Percentile p (0 < p < 100) of the n nonmissing values in work, which is
    reordered. by_formula selects egen's by-group arithmetic (see the
    definition); otherwise the rule of _pctile, which egen uses without by.
    n 0 gives missing.
*/
double cegen_percentile(double *work, size_t n, double p, int by_formula);

/* Median of the n nonmissing values in work (reordered); missing when n is 0 */
double cegen_median(double *work, size_t n);

/* ---------------------------------------------------------------------------
   Storage types and the store phase (cegen_store.c)
   --------------------------------------------------------------------------- */

/* Name to type; returns 0 and leaves *type alone when the name is unknown */
int cegen_type_parse(const char *name, cegen_type *type);
const char *cegen_type_name(cegen_type type);

/*
    The value a variable of the type holds after `generate type v = value':
    integer types truncate toward zero, float rounds to single precision, and
    a value outside the type's range becomes missing. Missing values keep
    their code.
*/
double cegen_type_round(double value, cegen_type type);

/*
    The type of a variable of the given type after replace has stored values
    whose nonmissing ones (any: there is one) lie in lo .. hi (fractional:
    some is not an integer). An integer type widens as replace widens it: to
    float (double from long) for a fraction, else to the narrowest integer
    type from its own up that holds the range, else to double. float and
    double stay. (cegen_rows.c)
*/
cegen_type cegen_chain_promote(cegen_type type, int any, int fractional, double lo, double hi);

/* Range and kind of the nonmissing values of a numeric result */
typedef struct
{
    size_t nonmissing;
    int all_integer; /* every nonmissing value is an integer */
    int all_float;   /* every nonmissing value is exactly a float */
    double min, max; /* of the nonmissing values (0 when there are none) */
} cegen_range;
void cegen_result_range(const double *values, size_t n, cegen_range *range);

/*
    Write a numeric result to variable var (one-based plugin position): every
    observation 1..nobs_total gets its row's value rounded to the type, or
    fill when it is not in obs_map. obs_map is strictly increasing. Counts
    the missing values written into *nmissing.
*/
stata_retcode cegen_store_numeric(const double *values, size_t n, const perm_idx_t *obs_map,
                                  int var, size_t nobs_total, cegen_type type, double fill,
                                  size_t *nmissing);

/* The same for a string result: "" outside the sample or for a NULL pointer */
stata_retcode cegen_store_strings(char *const *values, size_t n, const perm_idx_t *obs_map,
                                  int var, size_t nobs_total, size_t *nmissing);

/* Threads to use for a loop over n items (1 below the parallel threshold) */
int cegen_threads_for(size_t n);

#endif /* CEGEN_INTERNAL_H */

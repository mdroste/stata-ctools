/*
 * ccollapse_internal.h
 *
 * Types shared by the ccollapse modules:
 *   ccollapse_impl.c      arguments, the compute and store calls, the sample
 *   ccollapse_group.c     observations -> groups in the order of their keys
 *   ccollapse_stats.c     sums, moments, extremes, first and last values
 *   ccollapse_quantile.c  percentiles and the interquartile range
 *
 * Index spaces: rows are zero-based positions among the loaded (selected)
 * observations, in their order in the dataset; groups are numbered from zero
 * in the sort order of their by() values; sources and targets are zero based.
 */

#ifndef CCOLLAPSE_INTERNAL_H
#define CCOLLAPSE_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include "stplugin.h"
#include "ctools_data.h"

#define CCOLLAPSE_MODULE "ccollapse"

/* Statistics, as collapse names them. CC_STAT_PCTILE carries its percentile
 * (1..99) in cc_target.pct; median is p50. */
typedef enum
{
    CC_STAT_MEAN = 0,
    CC_STAT_SUM,
    CC_STAT_RAWSUM,
    CC_STAT_COUNT,
    CC_STAT_PERCENT,
    CC_STAT_SD,
    CC_STAT_SEMEAN,
    CC_STAT_SEBINOMIAL,
    CC_STAT_SEPOISSON,
    CC_STAT_MIN,
    CC_STAT_MAX,
    CC_STAT_FIRST,
    CC_STAT_LAST,
    CC_STAT_FIRSTNM,
    CC_STAT_LASTNM,
    CC_STAT_IQR,
    CC_STAT_PCTILE
} cc_stat;

typedef enum
{
    CC_W_NONE = 0,
    CC_W_FWEIGHT,
    CC_W_AWEIGHT,
    CC_W_IWEIGHT,
    CC_W_PWEIGHT
} cc_wtype;

typedef struct
{
    int source;   /* zero-based source variable */
    cc_stat stat; /* statistic */
    int pct;      /* percentile for CC_STAT_PCTILE, else 0 */
} cc_target;

/*
    The groups of the loaded rows.

    gid[i] is the group of row i; first[g] is the lowest row of group g, which
    supplies the by() values of the group. Groups are numbered in the sort
    order of their keys (numbers ascending with missing values last, strings
    by their bytes), which is the order of the collapsed dataset.

    Rows in group order are perm[starts[g] .. starts[g+1]), ascending within
    the group. perm == NULL means the rows already lie in group order
    (contiguous), so that position and row coincide. starts and perm are made
    by the grouping method that has them for free, and otherwise on demand by
    cc_groups_order().
*/
typedef struct
{
    size_t nobs;
    size_t ngroups;
    uint32_t *gid;      /* [nobs] */
    perm_idx_t *first;  /* [ngroups] */
    size_t *starts;     /* [ngroups + 1], or NULL until needed */
    perm_idx_t *perm;   /* [nobs], or NULL (identity when contiguous) */
    int contiguous;     /* rows are already in group order */
    const char *method; /* how the groups were found, for verbose */
    int method_code;    /* 1 one group, 2 ordered, 3 dense, 4 hash, 5 sort */
} cc_groups;

/*
    Find the groups of data->vars[0 .. nby). Every row belongs to a group.
    With nby == 0 all rows form one group.
*/
stata_retcode cc_group_rows(const stata_data *data, size_t nby, int threads, cc_groups *out);

/* Make out->starts (and out->perm unless contiguous) available. */
stata_retcode cc_groups_order(cc_groups *groups, int threads);

void cc_groups_free(cc_groups *groups);

/*
    One computation: the loaded columns, their groups and the targets.

    Numeric results go to values[t][g]; string results (first, last, firstnm,
    lastnm of a string source) are reported as rows: string_rows[t][g] is the
    row whose value group g takes. Exactly one of the two is non-NULL per
    target.
*/
typedef struct
{
    const stata_data *data;  /* loaded columns: by, sources, weight */
    size_t nby;
    size_t nsources;
    const double *weight;    /* [nobs], or NULL when unweighted */
    cc_wtype wtype;
    const cc_groups *groups;
    const cc_target *targets;
    size_t ntargets;
    double **values;         /* [ntargets] -> [ngroups], NULL for string targets */
    perm_idx_t **string_rows; /* [ntargets] -> [ngroups], NULL for numeric targets */
    int threads;
} cc_job;

/* Every statistic except the percentiles and iqr. */
stata_retcode cc_compute_stats(const cc_job *job);

/* The percentiles and iqr. Needs job->groups->starts (cc_groups_order). */
stata_retcode cc_compute_quantiles(const cc_job *job);

/* 1 for CC_STAT_PCTILE and CC_STAT_IQR */
static inline int cc_stat_is_quantile(cc_stat stat)
{
    return stat == CC_STAT_PCTILE || stat == CC_STAT_IQR;
}

/* A computed number as Stata stores it: missing when it is not a number, is
 * infinite, or lies in the range of the missing values. */
static inline double cc_storable(double value, double missval)
{
    return value < missval && value >= -1.7976931348623157e308 ? value : missval;
}

/* No row / no group */
#define CC_NONE UINT32_MAX

#endif /* CCOLLAPSE_INTERNAL_H */

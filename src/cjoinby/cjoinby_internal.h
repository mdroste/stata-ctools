/*
    cjoinby_internal.h
    Shared state and module interfaces of cjoinby

    Modules:
        cjoinby_table.c   staged observations (row-major, or columns when the
                          strings are too wide to stage) and their tile layout
        cjoinby_index.c   using-side key index: stable sort, key groups and a
                          direct or hashed lookup from key values to groups
        cjoinby_plan.c    master probe, match counts and the output layout
        cjoinby_write.c   fused, row-tiled writer of the joined observations
        cjoinby_impl.c    plugin calls, argument parsing, timing and cleanup

    Index spaces: plugin variable indices are one-based; rows, groups, drivers
    and output positions are zero-based; obs1 arguments are one-based Stata
    observations.
*/

#ifndef CJOINBY_INTERNAL_H
#define CJOINBY_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_rowstage.h"

/* No group, no master row, no using row */
#define CJ_NONE UINT32_MAX

/* unmatched(): which unmatched observations the output keeps */
#define CJ_KEEP_MASTER 1
#define CJ_KEEP_USING 2

/* _merge codes (joinby's): master only, using only, both */
#define CJ_MERGE_MASTER 1
#define CJ_MERGE_USING 2
#define CJ_MERGE_BOTH 3

/* ============================================================================
 * Staged observations
 *
 * Whole observations of a list of variables, held either row-major
 * (ctools_rowstage: one Stata record read per observation, one staged row
 * gathered per output row) or, when a string is too wide to stage, as loaded
 * columns. Either way a row is copied into a tile slot of the same layout:
 * offset[j] bytes into the slot, a double when width[j] is 0, otherwise a
 * NUL-terminated string in width[j] + 1 bytes.
 * ============================================================================ */

typedef struct
{
    int rowmajor; /* 1: rs holds the rows; 0: cols holds columns */
    ctools_rowstage rs;
    ctools_filtered_data cols;
    size_t nobs; /* staged observations */
    size_t nvars;
    size_t stride;   /* bytes per tile slot (multiple of 8) */
    int *var_idx;    /* [nvars] one-based plugin variable index */
    int *width;      /* [nvars] 0 numeric, else string width */
    size_t *offset;  /* [nvars] byte offset within a slot */
    int owns_layout; /* layout arrays are ours (column backend) */
} cj_table;

void cj_table_init(cj_table *t);
/* free now; release defers the large blocks to the release helper */
void cj_table_free(cj_table *t);
void cj_table_release(cj_table *t);

/* Stage observations obs1 .. obs1 + nobs - 1 of the listed variables (numeric
 * or str#; a strL fails with STATA_ERR_UNSUPPORTED_TYPE). nobs may be 0. */
stata_retcode cj_table_load(cj_table *t, const int *var_idx, size_t nvars, size_t obs1,
                            size_t nobs);

/* Write every staged observation back in order, starting at obs1. Allocates
 * nothing on the row backend, so it can undo a write that ran out of memory. */
stata_retcode cj_table_restore(const cj_table *t, size_t obs1);

/* Column col as a stata_variable for sorting and comparing: an owned copy on
 * the row backend, a borrowed view on the column backend. */
stata_retcode cj_table_column(const cj_table *t, size_t col, stata_variable *out);

/* Column-backend gather (out of line); see cj_table_gather */
void cj_table_gather_columns(const cj_table *t, size_t row, char *slot);

/* Copy staged row `row` into a tile slot of t->stride bytes */
static inline void cj_table_gather(const cj_table *t, size_t row, char *slot)
{
    if (t->rowmajor)
        memcpy(slot, t->rs.rows + row * t->stride, t->stride);
    else
        cj_table_gather_columns(t, row, slot);
}

/* A slot with no source row: missing numbers, empty strings */
void cj_table_blank(const cj_table *t, char *slot);

/* ============================================================================
 * Using-side key index
 * ============================================================================ */

/* The distinct key values, one entry per key group, in key order */
typedef struct
{
    int is_str;
    double *num; /* [ngroups] */
    char *str;   /* [ngroups * slot] NUL-terminated, zero padded */
    size_t slot; /* bytes per string entry */
} cj_keycol;

/* Direct lookup: every key is an integer (or a missing value) and the
 * product of the key ranges is small, so the group of a key tuple is one
 * table read. base/span/mult describe the mixed-radix position of each key. */
typedef struct
{
    uint32_t *table; /* [size] group, CJ_NONE when the tuple has no group */
    size_t size;
    double *base;     /* [nkeys] smallest finite value */
    double *top;      /* [nkeys] largest finite value */
    uint64_t *finite; /* [nkeys] number of finite codes (top - base + 1) */
    uint64_t *mult;   /* [nkeys] radix multiplier */
    int *has_missing; /* [nkeys] missing codes follow the finite ones */
} cj_direct;

/* Hashed lookup: open addressing, linear probing, one 8-byte slot per entry
 * (the group in the low word, 32 hash bits in the high word). */
typedef struct
{
    uint64_t *slots; /* [capacity], CJ_HASH_EMPTY when free */
    size_t mask;     /* capacity - 1 */
} cj_hash;

typedef struct
{
    int nkeys;
    size_t nobs;     /* using observations */
    stata_data keys; /* key columns; keys.sort_order is the stable key order */
    int *keycol;     /* [nkeys] 0 .. nkeys - 1 */
    size_t ngroups;
    uint32_t *gstart;  /* [ngroups + 1] first sorted position of each group */
    cj_keycol *keytab; /* [nkeys] */
    int use_direct;
    cj_direct direct;
    cj_hash hash;
} cj_index;

void cj_index_init(cj_index *ix);
void cj_index_free(cj_index *ix);
void cj_index_release(cj_index *ix);

/* Build the index over the first nkeys columns of utab. master_hint is the
 * master's observation count (sizes the direct table), 0 if unknown. */
stata_retcode cj_index_build(cj_index *ix, const cj_table *utab, int nkeys, size_t master_hint);

/* The group of every row of mkeys (its first nkeys columns, in key order):
 * gid[i] is the group or CJ_NONE; matched[g] is set to 1 for every group a row
 * found. Returns STATA_ERR_INVALID_INPUT when a key's type differs. */
stata_retcode cj_index_probe(const cj_index *ix, const stata_data *mkeys, uint32_t *gid,
                             unsigned char *matched);

/* Using rows of group g are order[gstart[g]] .. order[gstart[g + 1] - 1] */
static inline uint32_t cj_group_count(const cj_index *ix, uint32_t g)
{
    return ix->gstart[g + 1] - ix->gstart[g];
}

/* ============================================================================
 * Output plan
 *
 * The output is a sequence of drivers. A driver is a master observation
 * (paired with each using observation of its key group, or kept alone when it
 * has none and unmatched() keeps it) or a key group no master observation has
 * (unmatched(using|both): one output row per using observation). Driver d
 * owns output rows off[d] .. off[d + 1] - 1.
 * ============================================================================ */

typedef struct
{
    size_t nm;         /* master observations */
    size_t ndrv;       /* drivers: nm master rows, then the using-only groups */
    uint32_t *gid;     /* [ndrv] key group, CJ_NONE for an unmatched master row */
    uint32_t *off;     /* [ndrv + 1] first output row of each driver */
    uint32_t *row;     /* [ndrv] master row, CJ_NONE for a using-only group;
                          NULL when driver d < nm is master row d and the
                          using-only groups follow */
    size_t nout;       /* output observations */
    size_t n_master;   /* output rows from the master only */
    size_t n_using;    /* output rows from the using data only */
    size_t n_both;     /* matched output rows */
    size_t stage_from; /* output rows before this one are master rows that stay in place */
    int placed;        /* using-only groups were placed by the master's sort order */
} cj_plan;

void cj_plan_init(cj_plan *p);
void cj_plan_free(cj_plan *p);

/* Sort variables of the master, for placing using-only groups: plan variable
 * positions (zero-based columns of mkeys), keys first. */
typedef struct
{
    int nsort;
    const int *col;
} cj_sortspec;

/* Probe mkeys against the index and lay out the output. unmatched is a mask of
 * CJ_KEEP_*. With a sort spec the using-only groups are placed where a stable
 * sort on those variables would put them (master rows first among ties);
 * otherwise they follow the master rows, in key order. Returns a Stata return
 * code (920 when the output exceeds the observation limit). */
int cj_plan_build(cj_plan *p, const cj_index *ix, const stata_data *mkeys, size_t nm, int unmatched,
                  const cj_sortspec *sort);

static inline uint32_t cj_plan_master_row(const cj_plan *p, size_t d)
{
    return p->row ? p->row[d] : (d < p->nm ? (uint32_t)d : CJ_NONE);
}

/* ============================================================================
 * Writer
 * ============================================================================ */

typedef struct
{
    const cj_index *ix;
    const cj_plan *plan;
    const cj_table *utab; /* using rows: keys, then payload */
    const cj_table *mtab; /* master rows stage_from .. nm - 1, or NULL when none move */
    size_t stage_from;    /* first output row whose master variables are written */
    int nkeys;
    const int *key_mcol; /* [nkeys] zero-based column of each key in mtab's layout */
    const int *key_var;  /* [nkeys] one-based plugin index of each master key */
    int npay;
    const int *pay_mcol; /* [npay] zero-based mtab column of a shared variable, -1 for a new one */
    const int *pay_var;  /* [npay] one-based plugin index of a new variable, 0 for a shared one */
    int update, replace;
    int merge_var;          /* one-based, 0 none */
    int mobs_var, uobs_var; /* one-based, 0 none: source observation numbers */
    size_t limit;           /* write output rows [0, limit) only (fault hook), else nout */
} cj_writer;

stata_retcode cj_write_output(const cj_writer *w);

#endif /* CJOINBY_INTERNAL_H */

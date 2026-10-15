#ifndef CTOOLS_ROWSTAGE_H
#define CTOOLS_ROWSTAGE_H
#include "ctools_data.h"

/*
    Row-major staging for commands that move whole observations (csort, cmerge).

    Stata keeps a dataset as one record per observation, so reading or writing
    a variable for every observation walks the entire dataset once per
    variable, and a reordered write gathers each cell from a different place in
    its C column. Staging complete observations instead touches each Stata
    record once per direction and gathers one staged row per output row.
    Measured on 20 mixed variables (137-byte observations, 4 threads, Apple
    M4 Pro, October 2026): loads 2.1x and reordered stores 2.5x faster than
    the column transport at 10 million observations.

    Index spaces: var_idx[] holds one-based plugin variable indices; staged rows
    and source[] are zero-based; obs1 arguments are one-based Stata observations.
*/

/* A destination row with no staged source: numerics become missing, strings empty. */
#define CTOOLS_ROW_NONE PERM_IDX_MAX

typedef struct
{
    char *rows;     /* owned: nobs staged rows of stride bytes, plus a read guard */
    size_t stride;  /* bytes per staged row (multiple of 8) */
    size_t nobs;    /* staged rows */
    size_t nvars;   /* staged columns */
    int *var_idx;   /* owned [nvars]: one-based plugin variable index */
    int *width;     /* owned [nvars]: 0 = numeric (double), else fixed string width */
    size_t *offset; /* owned [nvars]: byte offset of the column within a row */
} ctools_rowstage;

void ctools_rowstage_init(ctools_rowstage *rs);
void ctools_rowstage_free(ctools_rowstage *rs);

/*
    1 when nobs observations of every listed variable can be staged: numeric,
    or a fixed-width string with a width hint from 1 to 2044 whose nobs slots
    fit in 2 GiB (the loader's bound for contiguous string columns; beyond it
    the column transport packs strings by their actual length). strL and
    unhinted strings cannot be staged.

    str_widths is indexed by plugin variable index minus one (SF_nvars()
    entries, as for ctools_data_load_ex); NULL reads the caller's _ctools_strw
    metadata.
*/
int ctools_rowstage_supported(const int *var_indices, size_t nvars, const int *str_widths,
                              size_t nobs);

/*
    Stage observations obs1 .. obs1+nobs-1 of the listed variables.

    Widths size the slots; every read is still bounded and a value longer than
    its hint fails the load.

    Returns STATA_ERR_UNSUPPORTED_TYPE before reading anything when a variable
    cannot be staged, so the caller can use the column transport instead.
    Nothing is left allocated on failure.
*/
stata_retcode ctools_rowstage_load(ctools_rowstage *rs, const int *var_indices, size_t nvars,
                                   const int *str_widths, size_t obs1, size_t nobs);

/*
    Copy staged column col into an owned, contiguous stata_variable (the layout
    ctools_data_load produces), for sorting or joining on it. Release with
    stata_variable_free or as part of a stata_data.
*/
stata_retcode ctools_rowstage_column(const ctools_rowstage *rs, size_t col, stata_variable *out);

/*
    Write every staged variable: observation obs1+d receives staged row
    source[d] for d in [0, n_out). source[d] may repeat rows or be
    CTOOLS_ROW_NONE; source == NULL writes the staged rows in order
    (n_out <= nobs) and allocates nothing, so it can restore the observations a
    failed write changed. Invalid arguments are rejected before the first write.
*/
stata_retcode ctools_rowstage_store(const ctools_rowstage *rs, const perm_idx_t *source,
                                    size_t n_out, size_t obs1);

#endif /* CTOOLS_ROWSTAGE_H */

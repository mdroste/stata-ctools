/*
    civreghdfe_impl.c
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    Implements 2SLS/IV regression with HDFE absorption. Reuses shared
    infrastructure for HDFE demeaning, singleton detection, and VCE computation.

    Algorithm:
    1. Read data: y, X_endog, X_exog, Z (instruments), FE vars, weights, cluster
    2. Detect and remove singletons
    3. Partial out FEs from all variables using CG solver
    4. First stage: Regress each endogenous var on [exogenous + instruments]
    5. Second stage: Regress y on [exogenous + predicted endogenous]
    6. Compute VCE (corrected for 2SLS)
    7. Store results to Stata
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "civreghdfe_impl.h"
#include "civreghdfe_result.h"
#include "civreghdfe_input.h"
#include "civreghdfe_boot.h"

#include "../ctools_config.h"
#include "../ctools_runtime.h"
#include "../ctools_ols.h"
#include "../ctools_matrix.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */
#include "civreghdfe_matrix.h"
#include "civreghdfe_estimate.h"

#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "../ctools_hdfe_utils.h"
#include "../hdfe/hdfe_sample.h"
#include "../creghdfe/creghdfe_hdfe.h"  /* heterogeneous slopes, shared with creghdfe */
#include "../ctools_unroll.h"
#include "../ctools_release.h"
#include "../ctools_threads.h"
#include "../hdfe/hdfe_solver.h"

static HDFE_State *iv_state = NULL;

/* vce(bootstrap): set by civreghdfe_main (civreghdfe_boot.c) */
ST_retcode (*civreghdfe_boot_retain_hook)(const civreghdfe_boot_source *src) = NULL;
void (*civreghdfe_boot_release_hook)(void) = NULL;
int civreghdfe_boot_retain_requested = 0;

/* Weighted sum of squares sum(w x^2); w may be NULL */
ST_double civ_wss(const ST_double *x, const ST_double *w, ST_int N)
{
    ST_double s = 0.0;
    if (w) { for (ST_int i = 0; i < N; i++) s += w[i] * x[i] * x[i]; }
    else   { for (ST_int i = 0; i < N; i++) s += x[i] * x[i]; }
    return s;
}

/* Weighted sum of squared deviations from the weighted mean; w may be NULL */
ST_double civ_wss_dev(const ST_double *x, const ST_double *w, ST_int N)
{
    ST_double sw = 0.0, m = 0.0, s = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_double wi = w ? w[i] : 1.0;
        sw += wi;
        m += wi * x[i];
    }
    if (sw <= 0.0) return 0.0;
    m /= sw;
    for (ST_int i = 0; i < N; i++) {
        ST_double d = x[i] - m;
        s += (w ? w[i] : 1.0) * d * d;
    }
    return s;
}

/* Clear keep[i] for every row whose value in v is missing (any missing code,
   or NaN). Same rows as SF_is_missing() for values read from Stata. */
static void civ_drop_missing(const ST_double *v, ST_int i0, ST_int i1, ST_double missval,
                             ST_int *keep)
{
    for (ST_int i = i0; i < i1; i++) {
        if (!(v[i] < missval)) keep[i] = 0;
    }
}

/* dst[j] = src[rows[j]], j < N: the retained rows of one column, in order;
   rows NULL: the first N rows as they are */
static void civ_gather(ST_double *dst, const ST_double *src, const ST_int *rows, ST_int N)
{
    if (!rows) {
        memcpy(dst, src, (size_t)N * sizeof(ST_double));
        return;
    }
    for (ST_int j = 0; j < N; j++) dst[j] = src[rows[j]];
}

/* The rows of src (N_total) where keep is nonzero, in order, into dst: row
   blocks in parallel, each starting where the kept rows before it end.
   Returns 0 or 920. */
static ST_retcode civ_gather_kept(ST_double *dst, const ST_double *src, const ST_int *keep,
                                  ST_int N_total)
{
    const int B = ctools_hdfe_row_blocks(N_total);
    ST_int *start = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    if (!start) return 920;
    int b;
    #pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++) {
        ST_int i0, i1, n = 0;
        ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
        for (ST_int ii = i0; ii < i1; ii++) n += keep[ii] != 0;
        start[b + 1] = n;
    }
    for (b = 0; b < B; b++) start[b + 1] += start[b];
    #pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++) {
        ST_int i0, i1, idx = start[b];
        ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
        for (ST_int ii = i0; ii < i1; ii++)
            if (keep[ii]) dst[idx++] = src[ii];
    }
    free(start);
    return 0;
}

/* The loaded rows that are valid and kept, ascending, into rows[0..N):
   kept holds the flags of the valid rows in order, or is NULL when all are
   kept. Row blocks in parallel, each starting where the valid rows and the
   kept rows before it end. Returns 0, 198 when the rows do not number N, or
   920. */
static ST_retcode civ_kept_rows(ST_int *rows, const ST_int *valid, const ST_int *kept,
                                ST_int N_total, ST_int N)
{
    const int B = ctools_hdfe_row_blocks(N_total);
    ST_int *valid_before = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    ST_int *kept_before = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    ST_retcode rc = 0;
    int b;
    if (!valid_before || !kept_before) {
        free(valid_before);
        free(kept_before);
        return 920;
    }
    #pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++) {
        ST_int i0, i1, n = 0;
        ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
        for (ST_int ii = i0; ii < i1; ii++) n += valid[ii] != 0;
        valid_before[b + 1] = n;
    }
    for (b = 0; b < B; b++) valid_before[b + 1] += valid_before[b];
    #pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++) {
        ST_int i0, i1, v = valid_before[b], n = 0;
        ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
        for (ST_int ii = i0; ii < i1; ii++) {
            if (!valid[ii]) continue;
            n += !kept || kept[v];
            v++;
        }
        kept_before[b + 1] = n;
    }
    for (b = 0; b < B; b++) kept_before[b + 1] += kept_before[b];
    if (kept_before[B] != N) {
        rc = 198;
    } else {
        #pragma omp parallel for schedule(static) if (B > 1)
        for (b = 0; b < B; b++) {
            ST_int i0, i1, v = valid_before[b], idx = kept_before[b];
            ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
            for (ST_int ii = i0; ii < i1; ii++) {
                if (!valid[ii]) continue;
                if (!kept || kept[v]) rows[idx++] = ii;
                v++;
            }
        }
    }
    free(valid_before);
    free(kept_before);
    return rc;
}

/*
    The columns the fixed effects are partialled out of.

    Each column needs, before the projection, its (weighted) sum of squares
    about zero and about its mean (civ_wss, civ_wss_dev) and its largest
    absolute value, the power of two below which it is divided by so that the
    solver works at unit scale; afterwards it is multiplied by that power
    again. Taken one after the other on a copy of the loaded column those are
    six passes over 80 MB per column at 10M rows, each a single chain of
    additions. Here the retained rows go from the loaded column straight to
    the column's place in the matrix that is partialled, and the passes are
    shared:

      civ_column_gather()  the rows, the largest absolute value and, without
                           weights, the sum and the sum of squares
      civ_column_sums()    with weights (normalized by then): their sum, the
                           weighted sum and the weighted sum of squares
      civ_column_center()  the sum of squared deviations, and the column
                           divided by its power of two
      civ_columns_restore() after the projection, every column times its
                           power of two again, on all threads

    With a single absorbed effect the projection is each value less the
    (weighted) mean of its level, and civ_column_absorb() takes the place of
    civ_column_center(), the projection and civ_columns_restore(): the level
    sums are added up from the scaled values in the pass that takes the
    squared deviations, and one more pass writes the column projected.

    Every sum adds the same values in row order as civ_wss(), civ_wss_dev()
    and project_and_subtract_fe() (hdfe_solver.c) do, a column on one thread,
    and every value is rounded by the same operations: the results are theirs
    to the last bit.
*/
typedef struct {
    const ST_double *src; /* the loaded column */
    const ST_int *rows;   /* the loaded row of each retained row; NULL: the first N as they are */
    ST_double *x;         /* N: the column's place in the matrix */
    ST_double *copy;      /* N or NULL: receives the partialled column as well */
    const ST_double *w;   /* N weights as the regression uses them, or NULL */
    ST_int N;
    ST_double largest;    /* max |x| */
    ST_double sw, sum;    /* sum of the weights (N without), (weighted) sum */
    ST_double ss_unc;     /* civ_wss() of the column */
    ST_double ss_pre;     /* civ_wss_dev() of the column */
    ST_double scale;      /* the power of two the column is divided by; 1: left as it is */
} civ_column;

static void civ_column_gather(civ_column *c, int take_sums)
{
    const ST_double *RESTRICT src = c->src;
    const ST_int *RESTRICT rows = c->rows;
    ST_double *RESTRICT x = c->x;
    const ST_int N = c->N;
    ST_double largest = 0.0, sum = 0.0, ss = 0.0;
    if (take_sums) {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = src[rows ? rows[i] : i];
            const ST_double a = fabs(v);
            x[i] = v;
            if (a > largest) largest = a;
            sum += v;
            ss += v * v;
        }
        c->sw = (ST_double)N;
        c->sum = sum;
        c->ss_unc = ss;
    } else {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = src[rows ? rows[i] : i];
            const ST_double a = fabs(v);
            x[i] = v;
            if (a > largest) largest = a;
        }
    }
    c->largest = largest;
}

static void civ_column_sums(civ_column *c)
{
    const ST_double *RESTRICT x = c->x, *RESTRICT w = c->w;
    const ST_int N = c->N;
    ST_double sw = 0.0, sum = 0.0, ss = 0.0;
    for (ST_int i = 0; i < N; i++) {
        const ST_double wx = w[i] * x[i];
        sw += w[i];
        sum += wx;
        ss += wx * x[i];
    }
    c->sw = sw;
    c->sum = sum;
    c->ss_unc = ss;
}

/* The power of two at or below the largest absolute value (1 for a column
   of zeros or with a value that is not finite) */
static void civ_column_set_scale(civ_column *c)
{
    int exponent = 0;
    c->scale = 1.0;
    if (c->largest > 0.0 && isfinite(c->largest)) {
        frexp(c->largest, &exponent);
        c->scale = ldexp(1.0, exponent - 1);
    }
}

/* x / scale. A power of two whose reciprocal is one too divides and
   multiplies to the same double; the others (a scale near the ends of the
   exponent range) are divided by. */
#define CIV_SCALED(v, scale, inv, by_inv) ((by_inv) ? (v) * (inv) : (v) / (scale))
static int civ_scale_inverse(ST_double scale, ST_double *inv)
{
    *inv = 1.0;
    if (!(scale >= 0x1p-1000 && scale <= 0x1p1000)) return 0;
    *inv = 1.0 / scale;
    return 1;
}

static void civ_column_center(civ_column *c)
{
    ST_double *RESTRICT x = c->x;
    const ST_double *RESTRICT w = c->w;
    const ST_int N = c->N;
    civ_column_set_scale(c);
    const ST_double scale = c->scale;
    ST_double inv;
    const int by_inv = civ_scale_inverse(scale, &inv);
    const int rescale = scale != 1.0;
    ST_double s = 0.0;
    if (c->sw <= 0.0) {
        c->ss_pre = 0.0;
        if (rescale) {
            for (ST_int i = 0; i < N; i++) x[i] = CIV_SCALED(x[i], scale, inv, by_inv);
        }
        return;
    }
    const ST_double m = c->sum / c->sw;
    if (w) {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = x[i], d = v - m;
            s += w[i] * d * d;
            if (rescale) x[i] = CIV_SCALED(v, scale, inv, by_inv);
        }
    } else {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = x[i], d = v - m;
            s += d * d;
            if (rescale) x[i] = CIV_SCALED(v, scale, inv, by_inv);
        }
    }
    c->ss_pre = s;
}

/* civ_column_center(), the projection on the one absorbed effect f (no
   slope) and the return to the column's units. sums: f->num_levels doubles
   of scratch. */
static void civ_column_absorb(civ_column *c, const FE_Factor *f, ST_double *RESTRICT sums)
{
    ST_double *RESTRICT x = c->x, *RESTRICT copy = c->copy;
    const ST_double *RESTRICT w = c->w;
    const ST_int *RESTRICT levels = f->levels;
    const ST_int N = c->N, num_levels = f->num_levels;
    /* the counts project_and_subtract_fe() divides by */
    const ST_double *RESTRICT inv_counts =
        (w != NULL && f->inv_weighted_counts != NULL) ? f->inv_weighted_counts : f->inv_counts;
    const ST_double *RESTRICT counts =
        (w != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;
    civ_column_set_scale(c);
    const ST_double scale = c->scale;
    ST_double inv;
    const int by_inv = civ_scale_inverse(scale, &inv);
    const int rescale = scale != 1.0;
    const int centered = !(c->sw <= 0.0);
    const ST_double m = centered ? c->sum / c->sw : 0.0;
    ST_double s = 0.0;
    memset(sums, 0, (size_t)num_levels * sizeof(ST_double));
    if (w) {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = x[i], d = v - m;
            const ST_double vs = rescale ? CIV_SCALED(v, scale, inv, by_inv) : v;
            s += w[i] * d * d;
            sums[levels[i] - 1] += vs * w[i];
        }
    } else {
        for (ST_int i = 0; i < N; i++) {
            const ST_double v = x[i], d = v - m;
            const ST_double vs = rescale ? CIV_SCALED(v, scale, inv, by_inv) : v;
            s += d * d;
            sums[levels[i] - 1] += vs;
        }
    }
    c->ss_pre = centered ? s : 0.0;
    if (inv_counts != NULL) {
        for (ST_int l = 0; l < num_levels; l++) sums[l] *= inv_counts[l];
    } else {
        for (ST_int l = 0; l < num_levels; l++) {
            if (counts[l] > 0) sums[l] /= counts[l];
        }
    }
    for (ST_int i = 0; i < N; i++) {
        const ST_double vs = rescale ? CIV_SCALED(x[i], scale, inv, by_inv) : x[i];
        ST_double r = vs - sums[levels[i] - 1];
        if (rescale) r *= scale;
        x[i] = r;
        if (copy) copy[i] = r;
    }
}

/* Rows per block of civ_columns_restore() */
#define CIV_RESTORE_BLOCK 16384

/* Every column times its power of two, and the columns that have one copied
   to their copy: blocks of rows on all threads */
static void civ_columns_restore(civ_column *cols, ST_int ncols, ST_int N)
{
    const ST_int nblk = N / CIV_RESTORE_BLOCK + (N % CIV_RESTORE_BLOCK != 0);
    ST_int blk;
    #pragma omp parallel for schedule(static) if ((size_t)N * (size_t)ncols > 200000)
    for (blk = 0; blk < nblk; blk++) {
        const ST_int i0 = blk * CIV_RESTORE_BLOCK;
        const ST_int i1 = (N - i0 < CIV_RESTORE_BLOCK) ? N : i0 + CIV_RESTORE_BLOCK;
        for (ST_int k = 0; k < ncols; k++) {
            ST_double *RESTRICT x = cols[k].x, *RESTRICT copy = cols[k].copy;
            const ST_double scale = cols[k].scale;
            if (scale != 1.0) {
                if (copy) {
                    for (ST_int i = i0; i < i1; i++) {
                        const ST_double v = x[i] * scale;
                        x[i] = v;
                        copy[i] = v;
                    }
                } else {
                    for (ST_int i = i0; i < i1; i++) x[i] *= scale;
                }
            } else if (copy) {
                memcpy(copy + i0, x + i0, (size_t)(i1 - i0) * sizeof(ST_double));
            }
        }
    }
}

/* A column's passes between the gather and the projection, as a task of the
   helper threads (the solver is prepared meanwhile) or of a loop */
static void *civ_column_center_run(void *arg)
{
    civ_column *c = (civ_column *)arg;
    if (c->w) civ_column_sums(c);
    civ_column_center(c);
    return NULL;
}

/* A loaded variable's memory back to the system, as a task of the helper
   threads (2 ms per 80 MB on Windows, and the next step need not wait) */
static void *civ_variable_free_run(void *arg)
{
    stata_variable_free((stata_variable *)arg);
    return NULL;
}

static void civ_release_prepared(void *prepared)
{
    partial_out_prepared_free((struct hdfe_coef_prepared *)prepared);
}

/* Values (rows times columns) from which the columns' passes run on the
   helper threads beside the solver's preparation */
#define CIV_AHEAD_MIN_VALUES ((size_t)1000000)

/* The slope description creghdfe's slope routines read: factor g is a
   slope when slope_idx[g] != 0, of the absorbed effect fe_absvar[g] */
static void civ_slope_input(const civreghdfe_input *in, creghdfe_input *ci)
{
    memset(ci, 0, sizeof(*ci));
    for (ST_int g = 0; g < 10; g++) {
        ci->slope_idx[g] = in->slope_flag[g];
        ci->fe_absvar[g] = in->fe_absvar[g];
    }
    ci->num_slopes = in->num_slopes;
    ci->dof_continuous = in->dof_continuous;
}

/* Redundant coefficients of each factor as reghdfe counts them (nested:
   NULL or G flags), and exact[g] where reghdfe reports the count as exact:
   nested effects, the first remaining intercept and, with mobility-group
   adjustments, the second (its "?" marks the others). With slopes the
   intercepts are counted among themselves and a slope constant within a
   level (zero without an intercept) is redundant under
   dofadjustments(continuous). Returns 0 or 920. */
static ST_retcode civ_redundant(const HDFE_State *S, const civreghdfe_input *in,
                                const ST_int *nested, ST_int dof_pairs, ST_int *redundant,
                                ST_int *exact)
{
    const ST_int G = S->G;
    if (in->num_slopes > 0 && G <= 10) {
        creghdfe_input ci;
        civ_slope_input(in, &ci);
        return creghdfe_redundant(S, &ci, nested, dof_pairs, redundant, exact) ? 920 : 0;
    }
    if (ctools_hdfe_redundant_intercepts(S->factors, G, S->N, nested, dof_pairs, redundant))
        return 920;
    ST_int first = -1, second = -1;
    for (ST_int g = 0; g < G; g++) {
        exact[g] = 0;
        if (nested && nested[g]) exact[g] = 1;
        else if (first < 0) first = g;
        else if (second < 0) second = g;
    }
    if (first >= 0) exact[first] = 1;
    if (dof_pairs && second >= 0) exact[second] = 1;
    return 0;
}

/*
    Read the 1 x n matrix of 0/1 column flags that the ado passes for
    partial(), orthog(), endogtest() and redundant() (one flag per expanded
    plugin column). Returns the number of flagged columns, or -1 if the
    matrix is missing or has the wrong shape.
*/
static ST_int read_flag_matrix(const char *name, ST_int n, ST_int *flags)
{
    char mname[64];
    snprintf(mname, sizeof(mname), "%s", name);
    if (n <= 0 || SF_row(mname) != 1 || SF_col(mname) != n) return -1;
    ST_int count = 0;
    for (ST_int k = 0; k < n; k++) {
        ST_double v = 0.0;
        if (SF_mat_el(mname, 1, k + 1, &v) != 0) return -1;
        flags[k] = (v != 0.0);
        count += flags[k];
    }
    return count;
}

/*
    Map the ado's flags for orthog()/endogtest()/redundant() (matrix `name`,
    one flag per column of the original list) onto the columns that survived
    the collinearity checks (dropped[j] != 0 marks a removed column; NULL =
    none removed). Writes the distinct 1-based positions among the surviving
    columns to out[0..n_flag-1]. A flagged column that was dropped cannot be
    tested (ivreg2 rejects it as well): its option code and position are
    saved for the ado's message and 198 is returned.
*/
static ST_retcode map_test_cols(const char *name, const char *opt, ST_int opt_code,
                                ST_int n_flag, ST_int n_orig, const ST_int *dropped,
                                ST_int *out)
{
    char msg[200];
    if (n_flag <= 0) return 0;
    ST_int *flags = (ST_int *)calloc((size_t)(n_orig > 0 ? n_orig : 1), sizeof(ST_int));
    if (!flags) {
        SF_error("civreghdfe: Memory allocation failed for test flags\n");
        return 920;
    }
    if (read_flag_matrix(name, n_orig, flags) != n_flag) {
        snprintf(msg, sizeof(msg), "civreghdfe: %s() column flags do not match the variable list\n", opt);
        SF_error(msg);
        free(flags);
        return 198;
    }
    ST_int pos = 0, n = 0;
    for (ST_int j = 0; j < n_orig; j++) {
        ST_int keep = dropped ? (dropped[j] == 0) : 1;
        pos += keep;
        if (!flags[j]) continue;
        if (!keep) {
            snprintf(msg, sizeof(msg),
                     "civreghdfe: a variable in %s() was dropped as collinear or absorbed by the fixed effects\n", opt);
            SF_error(msg);
            ST_retcode save_rc = ctools_scal_save("__civreghdfe_dropped_opt", (ST_double)opt_code);
            if (!save_rc) save_rc = ctools_scal_save("__civreghdfe_dropped_col", (ST_double)(j + 1));
            free(flags);
            return save_rc ? save_rc : 198;
        }
        out[n++] = pos;
    }
    free(flags);
    return 0;
}

/*
    FWL step for partial(). data (N x total_cols, column-major) is already
    FE-demeaned; is_partial[k] flags X_exog column exog_start + k, whose copy
    in Z is zexog_start + k. Every other column is residualized on the
    partial columns P by W-weighted least squares, (P'WP)^-1 P'Wx with the
    estimator's weights: ivreg2's s_partial after ivreghdfe's HDFE
    partial_out. Partial columns absorbed by the fixed effects (weighted SS
    <= absorb_tol times their pre-partialling SS in ss_exog) or collinear
    with earlier partial columns are left out of P and marked 2, as ivreg2
    drops collinear partial() variables; the rest keep flag 1. Messages
    only with report (the bootstrap replicates call it without).
*/
ST_retcode civ_fwl_partial_out(ST_double *data, ST_int N, ST_int total_cols,
                               ST_int exog_start, ST_int zexog_start, ST_int K_exog,
                               ST_int *is_partial, const ST_double *w,
                               const ST_double *ss_exog, ST_double absorb_tol, int report)
{
    ST_int n_cand = 0, n_keep = 0;
    for (ST_int k = 0; k < K_exog; k++) if (is_partial[k]) n_cand++;
    if (n_cand == 0) return 0;

    ST_int *cols = (ST_int *)malloc((size_t)n_cand * sizeof(ST_int));
    ST_int *coll = (ST_int *)calloc((size_t)n_cand, sizeof(ST_int));
    ST_double *P = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtP = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtP_inv = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtWX = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)total_cols, sizeof(ST_double));
    ST_double *coef = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)total_cols, sizeof(ST_double));
    ST_retcode rc = 0;

    if (!cols || !coll || !P || !PtP || !PtP_inv || !PtWX || !coef) {
        if (report) SF_error("civreghdfe: Memory allocation failed for FWL workspace\n");
        rc = 920;
        goto fwl_done;
    }

    /* P = partial columns that were not absorbed by the fixed effects */
    for (ST_int k = 0; k < K_exog; k++) {
        if (!is_partial[k]) continue;
        const ST_double *x = data + (size_t)(exog_start + k) * N;
        ST_double ss = civ_wss(x, w, N);
        if (ss < 1e-30 || ss <= absorb_tol * ss_exog[k]) {
            is_partial[k] = 2;
            continue;
        }
        memcpy(P + (size_t)n_keep * N, x, (size_t)N * sizeof(ST_double));
        cols[n_keep++] = k;
    }

    /* Drop partial columns collinear with earlier ones (column order, as
       ivreg2's _rmcollright), using the equilibrated P'WP */
    if (n_keep > 1) {
        ctools_matmul_atdb(P, P, w, N, n_keep, n_keep, PtP);
        if (detect_collinearity_scaled(PtP, n_keep, coll, 0) < 0) {
            if (report) SF_error("civreghdfe: Collinearity detection failed for partial()\n");
            rc = 920;
            goto fwl_done;
        }
        ST_int m = 0;
        for (ST_int a = 0; a < n_keep; a++) {
            if (coll[a]) {
                is_partial[cols[a]] = 2;
                continue;
            }
            if (m != a) {
                memcpy(P + (size_t)m * N, P + (size_t)a * N, (size_t)N * sizeof(ST_double));
                cols[m] = cols[a];
            }
            m++;
        }
        n_keep = m;
    }
    if (n_keep == 0) goto fwl_done;

    /* coef = (P'WP)^-1 P'W data for every column */
    ctools_matmul_atdb(P, P, w, N, n_keep, n_keep, PtP);
    memcpy(PtP_inv, PtP, (size_t)n_keep * n_keep * sizeof(ST_double));
    if (ctools_cholesky(PtP_inv, n_keep) != 0 ||
        ctools_invert_from_cholesky(PtP_inv, n_keep, PtP_inv) != 0) {
        if (report) SF_error("civreghdfe: cannot partial out the partial() variables (P'WP is singular)\n");
        rc = 498;
        goto fwl_done;
    }
    ctools_matmul_atdb(P, data, w, N, n_keep, total_cols, PtWX);
    ctools_matmul_ab(PtP_inv, PtWX, n_keep, n_keep, total_cols, coef);

    /* Residualize all columns except the partial columns and their Z copies */
    #pragma omp parallel for schedule(dynamic) if ((size_t)N * total_cols > 100000)
    for (ST_int c = 0; c < total_cols; c++) {
        if (c >= exog_start && c < exog_start + K_exog && is_partial[c - exog_start]) continue;
        if (c >= zexog_start && c < zexog_start + K_exog && is_partial[c - zexog_start]) continue;
        ST_double *x = data + (size_t)c * N;
        const ST_double *b = coef + (size_t)c * n_keep;
        for (ST_int i = 0; i < N; i++) {
            ST_double fitted = 0.0;
            for (ST_int p = 0; p < n_keep; p++) fitted += P[(size_t)p * N + i] * b[p];
            x[i] -= fitted;
        }
    }

fwl_done:
    free(cols); free(coll); free(P); free(PtP); free(PtP_inv); free(PtWX); free(coef);
    return rc;
}


/* A column is absorbed by the fixed effects (or by partial()) when its
   partialled (weighted) SS is at most CIV_FE_NOISE times its uncentered SS
   before partialling: only rounding noise is left, which the projection
   leaves at about 1e-30 of that SS. A column with genuine within variation
   above this is kept, as ivreghdfe keeps it (ivreghdfe has no rule of its
   own: it partials first and drops only what invsym() drops).
   CIV_FE_NOISE is defined in civreghdfe_impl.h. */

/* Columns dropped by Stata's invsym() of the symmetric K x K cross-product
   A (column-major), the rule behind ivreghdfe's ranks: a column is zeroed
   when its remaining pivot is not positive, below 1e-9 of its diagonal, or
   below 1e-19 (as ctools_invsym). The first n_ordered columns are swept in
   order, the rest on the largest remaining diagonal, as invsym() does.
   Columns with drop[j] already set are left out. Returns the number of
   columns newly dropped, or -1 if memory runs out. */
ST_int civ_invsym_drops(const ST_double *A, ST_int K, ST_int n_ordered, ST_int *drop)
{
    if (K <= 0) return 0;
    ST_double *W = (ST_double *)ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(ST_double));
    ST_int *done = (ST_int *)calloc((size_t)K, sizeof(ST_int));
    if (!W || !done) { free(W); free(done); return -1; }
    memcpy(W, A, (size_t)K * K * sizeof(ST_double));
    for (ST_int j = 0; j < K; j++) done[j] = drop[j] ? 1 : 0;
    ST_int n_new = 0;
    for (;;) {
        ST_int p = -1;
        for (ST_int j = 0; j < K; j++) {
            if (done[j]) continue;
            if (j < n_ordered) { p = j; break; }
            if (p < 0 || W[(size_t)j * K + j] > W[(size_t)p * K + p]) p = j;
        }
        if (p < 0) break;
        done[p] = 1;
        ST_double d = W[(size_t)p * K + p], d0 = A[(size_t)p * K + p];
        if (!(d > 0.0) || !(d0 > 0.0) || d < 1e-9 * d0 || d < 1e-19) {
            drop[p] = 1;
            n_new++;
            continue;
        }
        /* Sweep the remaining columns on p */
        for (ST_int i = 0; i < K; i++) {
            if (done[i]) continue;
            ST_double wip = W[(size_t)p * K + i] / d;
            for (ST_int j = i; j < K; j++) {
                if (done[j]) continue;
                W[(size_t)j * K + i] -= wip * W[(size_t)j * K + p];
                W[(size_t)i * K + j] = W[(size_t)j * K + i];
            }
        }
    }
    free(W);
    free(done);
    return n_new;
}

/*
    Full IV regression with HDFE.

    Variable layout from Stata:
    [y, X_endog_1, ..., X_endog_Ke, X_exog_1, ..., X_exog_Kx, Z_1, ..., Z_Kz, FE_1, ..., FE_G, cluster, weight,
     tsset panel, tsset time]

    Scalars from Stata:
    - __civreghdfe_K_endog: number of endogenous regressors
    - __civreghdfe_K_exog: number of exogenous regressors
    - __civreghdfe_K_iv: number of instruments (including exogenous)
    - __civreghdfe_G: number of FE factors
    - __civreghdfe_has_cluster, __civreghdfe_has_weights, __civreghdfe_weight_type
    - __civreghdfe_vce_type: 0=unadjusted, 1=robust, 2=cluster
    - __civreghdfe_maxiter, __civreghdfe_tolerance, __civreghdfe_verbose
*/
/*
    ivreg2 options read here for the estimation (civreghdfe_extra):
    smatrix()/wmatrix() arrive as __civreghdfe_S0/__civreghdfe_W0 over the
    instruments left after partial() (Z's order: included, then excluded)
    and are cut to the instruments kept; bw(auto) needs the row sums of
    the unpartialled instruments (no constant, no partialled regressor);
    orthog() on included instruments arrives as __civreghdfe_orthogx_flags
    over the exogenous regressors left after partial().
*/
typedef struct {
    civreghdfe_extra ex;
    ST_double *S0, *W0, *abw_z;
    ST_int *ox;
    ST_int n_ox;
} civ_extra_buf;

static void civ_extra_free(civ_extra_buf *b)
{
    free(b->S0); free(b->W0); free(b->abw_z); free(b->ox);
    memset(b, 0, sizeof(*b));
}

static ST_retcode civ_read_user_matrix(const char *name, ST_int Lp, const ST_int *keep,
                                       ST_int K_iv, ST_double **out)
{
    char mname[64];
    snprintf(mname, sizeof(mname), "%s", name);
    *out = NULL;
    if (SF_row(mname) != Lp || SF_col(mname) != Lp) return 0;
    ST_double *M = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
    if (!M) return 920;
    ST_int c = 0;
    for (ST_int j = 0; j < Lp; j++) {
        if (!keep[j]) continue;
        ST_int r = 0;
        for (ST_int i = 0; i < Lp; i++) {
            if (!keep[i]) continue;
            ST_double v = 0.0;
            if (SF_mat_el(mname, i + 1, j + 1, &v) != 0 || !(v < SV_missval)) { free(M); return 504; }
            M[(size_t)c * K_iv + r] = v;
            r++;
        }
        c++;
    }
    *out = M;
    return 0;
}

static ST_retcode civ_build_extra(civ_extra_buf *b, ST_int N, ST_int K_exog_in, ST_int K_excl_in,
                                  const ST_double *Z_raw, ST_int G, ST_int K_exog_post,
                                  const ST_int *is_collinear_x, ST_int K_iv,
                                  ST_int kernel_type, ST_double bw, const ST_double *abw_fe)
{
    ST_double dv = 0;
    ST_retcode rc = 0;
    memset(b, 0, sizeof(*b));
    if (SF_scal_use("__civreghdfe_ffirst", &dv) == 0 && dv != 0) b->ex.ffirst = 1;
    if (SF_scal_use("__civreghdfe_psd", &dv) == 0 && (dv == 1 || dv == 2)) b->ex.psd = (ST_int)dv;
    if (SF_scal_use("__civreghdfe_store_sw", &dv) == 0 && dv != 0) b->ex.store_sw = 1;

    /* partial() flags of the raw exogenous columns */
    ST_int *part = (ST_int *)calloc((size_t)(K_exog_in > 0 ? K_exog_in : 1), sizeof(ST_int));
    if (!part) return 920;
    char mname[64];
    snprintf(mname, sizeof(mname), "__civreghdfe_partial_flags");
    if (SF_scal_use("__civreghdfe_n_partial", &dv) == 0 && dv > 0 &&
        SF_row(mname) == 1 && SF_col(mname) == K_exog_in) {
        for (ST_int k = 0; k < K_exog_in; k++) {
            ST_double v = 0;
            if (SF_mat_el(mname, 1, k + 1, &v) == 0) part[k] = (v != 0.0);
        }
    }

    /* kept instruments, over the columns left after partial() */
    ST_int Lp = K_exog_post + K_excl_in;
    ST_int *keep = (ST_int *)calloc((size_t)(Lp > 0 ? Lp : 1), sizeof(ST_int));
    if (!keep) { free(part); return 920; }
    for (ST_int k = 0; k < K_exog_post; k++) keep[k] = !is_collinear_x[k];
    for (ST_int j = 0; j < K_excl_in; j++) {
        char sn[64];
        snprintf(sn, sizeof(sn), "__civreghdfe_zcollin_%d", (int)(j + 1));
        ST_double z = 0;
        keep[K_exog_post + j] = !(SF_scal_use(sn, &z) == 0 && z != 0);
    }
    rc = civ_read_user_matrix("__civreghdfe_S0", Lp, keep, K_iv, &b->S0);
    if (!rc) rc = civ_read_user_matrix("__civreghdfe_W0", Lp, keep, K_iv, &b->W0);
    free(keep);

    /* bw(auto) */
    if (!rc && kernel_type > 0 && bw < 0 && (Z_raw || abw_fe)) {
        ST_double has_cons = 0, span = 0;
        SF_scal_use("__civreghdfe_has_cons", &has_cons);
        SF_scal_use("__civreghdfe_tspan", &span);
        ST_int cons_col = (G == 0 && has_cons != 0) ? K_exog_in - 1 : -1;
        b->abw_z = (ST_double *)calloc((size_t)(N > 0 ? N : 1), sizeof(ST_double));
        if (!b->abw_z) rc = 920;
        for (ST_int k = 0; !rc && !abw_fe && k < K_exog_in + K_excl_in; k++) {
            if (k < K_exog_in && (k == cons_col || part[k])) continue;
            const ST_double *z = Z_raw + (size_t)k * N;
            for (ST_int i = 0; i < N; i++) b->abw_z[i] += z[i];
        }
        if (abw_fe) memcpy(b->abw_z, abw_fe, (size_t)N * sizeof(ST_double));
        b->ex.abw_tobs = span + 1.0;
    }
    free(part);

    /* orthog() on included instruments */
    if (!rc && SF_scal_use("__civreghdfe_n_orthogx", &dv) == 0 && dv > 0) {
        ST_int n = (ST_int)dv;
        b->ox = (ST_int *)malloc((size_t)n * sizeof(ST_int));
        if (!b->ox) rc = 920;
        else {
            rc = map_test_cols("__civreghdfe_orthogx_flags", "orthog", 4, n, K_exog_post,
                               is_collinear_x, b->ox);
            b->n_ox = n;
        }
    }
    b->ex.S0 = b->S0;
    b->ex.W0 = b->W0;
    b->ex.abw_z = b->abw_z;
    if (rc == 504) SF_error("civreghdfe: smatrix()/wmatrix() has missing values\n");
    return rc;
}

/* The fit's HDFE state, released by the helper thread (ctools_release.h) */
static void civ_release_state(void *state)
{
    ctools_hdfe_state_cleanup((HDFE_State *)state);
    free(state);
}

static ST_retcode do_iv_regression(void)
{
    /* Heap state released at iv_cleanup; every error path leaves through it */
    ST_retcode iv_rc = 0;
    ST_int *cluster_ids = NULL, *cluster2_ids = NULL, *valid_mask = NULL, *singleton_mask = NULL;
    ST_int **fe_levels_c = NULL;
    /* What the FE remap and singleton removal learn about the factors, for
     * the steps after them (hdfe_sample.h) */
    HDFE_FactorOrder fe_order[10];
    HDFE_Singletons fe_sing;
    memset(fe_order, 0, sizeof(fe_order));
    memset(&fe_sing, 0, sizeof(fe_sing));
    ST_int *level_max = NULL, *is_partial = NULL, *is_collinear_x = NULL, *test_idx = NULL, *ts_int = NULL;
    ST_int *fe_redundant = NULL, *fe_nested = NULL, *fe_exact = NULL;
    ST_double *y_c = NULL, *X_endog_c = NULL, *X_exog_c = NULL, *Z_c = NULL, *weights_c = NULL;
    ST_double *slope_c[10] = {0};  /* compacted slope variables (slope factors only) */
    ST_int *cluster_ids_c = NULL, *cluster2_ids_c = NULL;
    ST_double *all_data = NULL, *beta = NULL, *V = NULL, *first_stage_F = NULL, *abw_fe = NULL;
    civ_column *part_cols = NULL;             /* the partialled columns of all_data */
    struct hdfe_coef_prepared *prepared = NULL; /* the projection's part done ahead of the columns */
    ctools_persistent_pool *free_pool = NULL; /* helper threads returning the loaded variables */
    perm_idx_t *est_obs = NULL;
    HDFE_State *state = NULL;
    /* vce(bootstrap): what the retained model needs beyond the fit's arrays */
    const int bs_retain = civreghdfe_boot_retain_requested && civreghdfe_boot_retain_hook;
    ST_int *bs_partial = NULL, *bs_collin_z = NULL;
    ST_double *bs_tilde = NULL, *bs_b = NULL;
    ST_int bs_endo_in_z = 0;

    /* Start total timing */
    double t_total_start = ctools_timer_seconds();
    double t_load = 0, t_extract = 0, t_singleton = 0, t_hdfe_setup = 0, t_fwl = 0, t_partial_out = 0;
    double t_postproc = 0, t_estimate = 0, t_store = 0;

    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N_total = in2 - in1 + 1;

    civreghdfe_input input;
    ST_retcode input_rc = civreghdfe_read_input(&input);
    if (input_rc) return input_rc;
    ST_int K_endog = input.K_endog;
    ST_int K_exog = input.K_exog;
    ST_int K_iv = input.K_iv;
    ST_int G = input.G;
    ST_int has_cluster = input.has_cluster;
    ST_int has_cluster2 = input.has_cluster2;
    ST_int has_weights = input.has_weights;
    ST_int weight_type = input.weight_type;
    ST_int vce_type = input.vce_type;
    ST_int maxiter = input.maxiter;
    ST_int verbose = input.verbose;
    ST_int est_method = input.est_method;
    ST_int kernel_type = input.kernel_type;
    ST_double bw = input.bw;
    ST_int kiefer = input.kiefer;
    ST_int ts_mode = input.ts_mode;
    ST_int ts_cols = input.ts_cols;
    ST_int dofminus = input.dofminus;
    ST_int sdofminus_opt = input.sdofminus_opt;
    ST_int nopartialsmall = input.nopartialsmall;
    ST_int center = input.center;
    ST_double tolerance = input.tolerance;
    ST_double kclass_user = input.kclass_user;
    ST_double fuller_alpha = input.fuller_alpha;
    ST_double tdelta = input.tdelta;

    /* The IV code indexes its K x K workspaces in int: j * K + i overflows
       once K exceeds 46340 (cqreg has the same cap) */
    if (K_iv > 46340 || K_exog + K_endog > 46340) {
        SF_error("civreghdfe: at most 46340 regressors and instruments are supported\n");
        return 198;
    }

    /* Calculate variable positions */
    /* Layout: [y, X_endog (Ke), X_exog (Kx), excluded instruments (Kz - Kx), FE (G),
                cluster?, cluster2?, weight?, tsset panel?, tsset time?] */
    ST_int var_y_idx = 0;  /* Index in loaded_data.vars[] */
    ST_int var_endog_start_idx = 1;
    ST_int var_exog_start_idx = var_endog_start_idx + K_endog;
    /* K_iv counts the exogenous regressors (the first columns of Z) and the
       excluded instruments; only the latter are sent a second time */
    const ST_int K_excl_in = K_iv - K_exog;
    ST_int var_excl_start_idx = var_exog_start_idx + K_exog;
    ST_int var_fe_start_idx = var_excl_start_idx + K_excl_in;
    ST_int var_cluster_idx = has_cluster ? (var_fe_start_idx + G) : -1;
    ST_int var_cluster2_idx = has_cluster2 ? (var_fe_start_idx + G + (has_cluster ? 1 : 0)) : -1;
    ST_int var_weight_idx = has_weights ? (var_fe_start_idx + G + (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0)) : -1;
    /* tsset columns (panel, time) follow the weight */
    ST_int var_ts_start_idx = var_fe_start_idx + G + (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0) + (has_weights ? 1 : 0);
    /* Heterogeneous slopes: one variable per slope factor, after the tsset columns */
    const ST_int num_slopes = input.num_slopes;
    ST_int var_slope_start_idx = var_ts_start_idx + ts_cols;

    ST_int K_total = K_exog + K_endog;  /* Total regressors */

    /* Start data load timing */
    double t_load_start = ctools_timer_seconds();

    /* ================================================================
     * PARALLEL DATA LOADING using ctools_data_load
     * This handles if/in filtering at load time, loading only filtered observations.
     * ================================================================ */

    /* Calculate total variables to load */
    if (K_excl_in < 0) {
        SF_error("civreghdfe: fewer instruments than exogenous regressors\n");
        return 198;
    }
    ST_int total_vars = 1 + K_endog + K_exog + K_excl_in + G +
                        (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0) + (has_weights ? 1 : 0) + ts_cols +
                        num_slopes;

    /* Build variable indices array (1-based for Stata) */
    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("civreghdfe: Memory allocation failed\n");
        return 920;
    }
    for (ST_int i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;  /* 1-based Stata variable indices */
    }

    /* Load all variables in parallel with if/in filtering */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER | (input.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, load_flags);
    free(var_indices);

    /* Capture pure SPI data load time */
    t_load = ctools_timer_seconds() - t_load_start;

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("civreghdfe: Parallel data load failed\n");
        return ctools_stata_rc(load_rc);
    }

    /* Start extraction timer */
    double t_extract_start = ctools_timer_seconds();

    /* Get filtered observation count - data is already filtered */
    N_total = (ST_int)filtered.data.nobs;

    if (N_total == 0) {
        SF_error("civreghdfe: No observations selected\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* Allocate cluster ID arrays (needed for string conversion + sentinel marking).
     * All other data (y, X, Z, FE, weights) is read directly from filtered.data
     * during missing check and compaction, avoiding redundant N_total-sized copies. */
    cluster_ids = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;
    cluster2_ids = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;

    if ((has_cluster && !cluster_ids) || (has_cluster2 && !cluster2_ids)) {
        SF_error("civreghdfe: Memory allocation failed\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Extract cluster IDs (needs early processing for string variables) */
    if (has_cluster) {
        if (filtered.data.vars[var_cluster_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster_idx].data.str,
                                               N_total, cluster_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster variable to group IDs\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
        } else {
            double *src = filtered.data.vars[var_cluster_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster_ids, &n_groups) != 0) {
                iv_rc = 920;
                goto iv_cleanup;
            }
        }
    }
    if (has_cluster2) {
        if (filtered.data.vars[var_cluster2_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster2_idx].data.str,
                                               N_total, cluster2_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster2 variable to group IDs\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
        } else {
            double *src = filtered.data.vars[var_cluster2_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster2_ids, &n_groups) != 0) {
                iv_rc = 920;
                goto iv_cleanup;
            }
        }
    }

    /* End extraction, start missing value check timing.
     * Data arrays (y, X, Z, FE, weights) are NOT copied — we read directly
     * from filtered.data during missing value check and compaction. */
    t_extract = ctools_timer_seconds() - t_extract_start;
    double t_missing_start = ctools_timer_seconds();

    /* Drop observations with missing values. Columns are scanned one at a
     * time over blocks of rows (sequential reads; blocks in parallel),
     * directly from filtered.data. */
    valid_mask = (ST_int *)malloc((size_t)N_total * sizeof(ST_int));
    if (!valid_mask) {
        SF_error("civreghdfe: Memory allocation failed for valid_mask\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    ST_int N_valid = 0;
    {
        const ST_double missval = SV_missval;
        /* y, X_endog, X_exog, excluded instruments and FE variables are the
           first var_fe_start_idx + G columns */
        const ST_int n_lead = var_fe_start_idx + G;
        const ST_int block = 16384;
        /* N_total + block - 1 would overflow int near INT_MAX rows */
        const ST_int num_blocks = N_total / block + (N_total % block != 0);
        ST_int blk;
        #pragma omp parallel for schedule(static) reduction(+ : N_valid) if (N_total > 100000)
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_int i0 = blk * block;
            const ST_int i1 = (N_total - i0 < block) ? N_total : i0 + block;
            for (ST_int ii = i0; ii < i1; ii++) valid_mask[ii] = 1;
            for (ST_int c = 0; c < n_lead; c++)
                civ_drop_missing(filtered.data.vars[c].data.dbl, i0, i1, missval, valid_mask);
            for (ST_int k = 0; k < ts_cols; k++)
                civ_drop_missing(filtered.data.vars[var_ts_start_idx + k].data.dbl, i0, i1, missval, valid_mask);
            for (ST_int k = 0; k < num_slopes; k++)
                civ_drop_missing(filtered.data.vars[var_slope_start_idx + k].data.dbl, i0, i1, missval, valid_mask);
            if (has_weights) {
                const ST_double *w = filtered.data.vars[var_weight_idx].data.dbl;
                for (ST_int ii = i0; ii < i1; ii++) {
                    if (!(w[ii] < missval) || w[ii] <= 0) valid_mask[ii] = 0;
                }
            }
            /* missing cluster values were marked -1 during extraction */
            if (has_cluster) {
                for (ST_int ii = i0; ii < i1; ii++) if (cluster_ids[ii] < 0) valid_mask[ii] = 0;
            }
            if (has_cluster2) {
                for (ST_int ii = i0; ii < i1; ii++) if (cluster2_ids[ii] < 0) valid_mask[ii] = 0;
            }
            for (ST_int ii = i0; ii < i1; ii++) N_valid += valid_mask[ii];
        }
    }

    /* Factor-variable indicators (the ado's _ctools_fvnative): with
     * __civreghdfe_fvspec, loaded column c of [X_endog, X_exog, excluded
     * instruments] (c = 1 ...) holds a factor variable when entry c of
     * __civreghdfe_fvlev is a level (>= 0). Now that its missing values have
     * been screened it becomes 1 where it equals that level, 0 elsewhere:
     * the doubles fvrevar's indicator holds, for every later use. */
    {
        ST_double use_fv = 0;
        const ST_int n_x = var_fe_start_idx - 1;
        if (SF_scal_use("__civreghdfe_fvspec", &use_fv) == 0 && use_fv == 1 && n_x > 0) {
            for (ST_int c = 1; c <= n_x; c++) {
                ST_double level = -1;
                if (SF_mat_el("__civreghdfe_fvlev", 1, c, &level) != 0 ||
                    !(level >= -1 && level == floor(level))) {
                    SF_error("civreghdfe: invalid factor-variable specification\n");
                    iv_rc = 198;
                    goto iv_cleanup;
                }
                if (level < 0) continue;
                ST_double *col = filtered.data.vars[c].data.dbl;
                ST_int r;
                #pragma omp parallel for schedule(static) if (N_total > 100000)
                for (r = 0; r < N_total; r++) col[r] = col[r] == level ? 1.0 : 0.0;
            }
            if ((iv_rc = ctools_scal_save("__civreghdfe_fvdone", 1.0))) goto iv_cleanup;
        }
    }

    if (N_valid <= K_total || N_valid < K_iv) {
        SF_error("civreghdfe: Insufficient observations\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* Compact FE levels from raw doubles with O(N) counting-based remap.
     * This replaces the old approach of: remap_values_sorted (O(N log N)) +
     * separate integer compaction. Now we compact and remap in one step. */
    fe_levels_c = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    int fe_c_alloc_failed = 0;
    if (fe_levels_c) {
        for (ST_int g = 0; g < G; g++) {
            fe_levels_c[g] = (ST_int *)malloc((size_t)N_valid * sizeof(ST_int));
            if (!fe_levels_c[g]) fe_c_alloc_failed = 1;
        }
    } else {
        fe_c_alloc_failed = 1;
    }
    if (fe_c_alloc_failed) {
        SF_error("civreghdfe: Memory allocation failed for FE compaction\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    {
        /* The valid rows of each FE variable (compacted in parallel row
         * blocks when some rows are missing) to levels 1..L in order of
         * first appearance, remapped in parallel, with the rows alone in a
         * level noted for singleton removal. level_max[g] keeps L, a bound
         * on the levels left then. */
        const int all_valid = N_valid == N_total;
        const int B = ctools_hdfe_row_blocks(N_total);
        double *fe_compact = all_valid ? NULL : (double *)malloc((size_t)N_valid * sizeof(double));
        ST_int *block_start = all_valid ? NULL : (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
        level_max = (ST_int *)calloc((size_t)G + 1, sizeof(ST_int));
        if (!level_max || (!all_valid && (!fe_compact || !block_start))) {
            free(fe_compact);
            free(block_start);
            SF_error("civreghdfe: Memory allocation failed for FE compaction\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
        if (!all_valid) {
            int b;
            #pragma omp parallel for schedule(static) if (B > 1)
            for (b = 0; b < B; b++) {
                ST_int i0, i1, n = 0;
                ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
                for (ST_int ii = i0; ii < i1; ii++) n += valid_mask[ii] != 0;
                block_start[b + 1] = n;
            }
            for (b = 0; b < B; b++) block_start[b + 1] += block_start[b];
        }

        for (ST_int g = 0; g < G; g++) {
            const double *fe_src = filtered.data.vars[var_fe_start_idx + g].data.dbl;
            if (!all_valid) {
                int b;
                #pragma omp parallel for schedule(static) if (B > 1)
                for (b = 0; b < B; b++) {
                    ST_int i0, i1, idx = block_start[b];
                    ctools_hdfe_block_range(N_total, B, b, &i0, &i1);
                    for (ST_int ii = i0; ii < i1; ii++)
                        if (valid_mask[ii]) fe_compact[idx++] = fe_src[ii];
                }
            }
            const int remap_rc = hdfe_remap_levels_ordered(all_valid ? fe_src : fe_compact, N_valid,
                                                           fe_levels_c[g], &level_max[g],
                                                           &fe_order[g].rows);
            if (remap_rc < 0) {
                SF_error("civreghdfe: FE level remapping failed\n");
                free(fe_compact);
                free(block_start);
                iv_rc = 920;
                goto iv_cleanup;
            }
            fe_order[g].ordered = remap_rc == 0;
        }
        free(fe_compact);
        free(block_start);
    }

    /* End missing value check/compact timing, start singleton timing */
    double t_missing = ctools_timer_seconds() - t_missing_start;
    double t_singleton_start = ctools_timer_seconds();

    /* Step 2: Singleton detection on compacted FE levels, from the rows the
     * remap found alone in a level. fe_sing.keep holds the flags of the
     * valid rows (NULL when none is dropped) until the factors are built */
    /* ivreghdfe selects the sample with reghdfe: fweights drop only levels
     * whose total fweight is 1; aweights and pweights use the unweighted rule. */
    /* the fweights of the valid rows: the loaded column itself when every
       row is valid */
    ST_double *fw_c = NULL;
    const ST_double *fw_rows = NULL;
    if (has_weights && weight_type == 2) {
        const ST_double *fw_src = filtered.data.vars[var_weight_idx].data.dbl;
        if (N_valid == N_total) fw_rows = fw_src;
        else {
            fw_c = (ST_double *)malloc((size_t)N_valid * sizeof(ST_double));
            if (fw_c && civ_gather_kept(fw_c, fw_src, valid_mask, N_total) != 0) {
                free(fw_c);
                fw_c = NULL;
            }
            fw_rows = fw_c;
        }
    }
    /* keepsingletons: reghdfe's option, which ivreghdfe passes on */
    ST_double keep_singletons = 0;
    SF_scal_use("__civreghdfe_keepsingletons", &keep_singletons);
    ST_int num_singletons_total;
    if (keep_singletons == 1) {
        num_singletons_total = 0;
    } else if ((has_weights && weight_type == 2 && !fw_rows) ||
               hdfe_drop_singletons(&fe_sing, fe_levels_c, level_max, fe_order, G, N_valid, fw_rows,
                                    verbose)) {
        num_singletons_total = -1;
    } else {
        num_singletons_total = fe_sing.num_dropped;
    }
    free(fw_c);
    if (num_singletons_total >= 0 &&
        (iv_rc = ctools_scal_save("__civreghdfe_num_singletons", (ST_double)num_singletons_total)))
        goto iv_cleanup;
    if (num_singletons_total < 0) {
        SF_error("civreghdfe: Memory allocation failed in singleton detection\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Step 3: Build combined index mapping original obs → final position.
     * Combine valid_mask and singleton_mask into one compaction pass. */
    ST_int N = N_valid - num_singletons_total;
    /* The estimation rows: with fweights __civreghdfe_N is their total
     * weight, and the wrapper's e(N_full) counts rows */
    if ((iv_rc = ctools_scal_save("__civreghdfe_N_rows", (ST_double)N))) goto iv_cleanup;

    /* The matrix that the fixed effects are partialled out of, N x
       total_cols column-major: Z = [copy of X_exog | excluded instruments],
       then y, X_endog and X_exog. The retained rows of every column go
       straight to its place in it; Z's copy of X_exog is filled from the
       partialled X_exog, which gives exactly what partialling it separately
       would. Two spare rows (total_cols values each after the N x total_cols
       block) hold each column's pre-partialling sums of squares, centred and
       uncentred; they share all_data's lifetime. */
    const ST_int total_cols = 1 + K_endog + K_exog + K_iv;
    const ST_int col_excl = K_exog;
    const ST_int col_y = K_iv;
    const ST_int col_endog = K_iv + 1;
    const ST_int col_exog = K_iv + 1 + K_endog;
    /* columns partialled: every one from the excluded instruments on */
    const ST_int n_part = total_cols - col_excl;
    all_data = (ST_double *)ctools_safe_malloc3((size_t)N + 2, (size_t)total_cols, sizeof(ST_double));
    part_cols = (civ_column *)calloc((size_t)n_part, sizeof(civ_column));
    /* The unpartialled columns are read again by the bootstrap's refits
       and, without absorbed effects, by bw(auto): only then are they kept
       apart, as y_c, X_endog_c, X_exog_c and Z_c = [X_exog | excluded
       instruments]. Z_c's storage ends with the ts_cols tsset columns
       (panel, time), which share its lifetime; without the copies it holds
       those alone. */
    const int keep_raw = bs_retain || (G == 0 && kernel_type > 0 && bw < 0);
    if (keep_raw) {
        y_c = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
        X_endog_c = (K_endog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_endog, sizeof(ST_double)) : NULL;
        X_exog_c = (K_exog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_exog, sizeof(ST_double)) : NULL;
    }
    const ST_int z_cols = (keep_raw ? K_iv : 0) + ts_cols;
    Z_c = z_cols > 0 ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)z_cols, sizeof(ST_double)) : NULL;
    ST_double *ts_c = Z_c ? Z_c + (size_t)N * (keep_raw ? K_iv : 0) : NULL;
    weights_c = has_weights ? (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double)) : NULL;
    cluster_ids_c = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;
    cluster2_ids_c = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;

    if (!all_data || !part_cols || (z_cols > 0 && !Z_c) ||
        (keep_raw && (!y_c || (K_endog > 0 && !X_endog_c) || (K_exog > 0 && !X_exog_c))) ||
        (has_weights && !weights_c) || (has_cluster && !cluster_ids_c) ||
        (has_cluster2 && !cluster2_ids_c)) {
        SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Compaction: the retained rows (valid and not singleton) of every
     * column, read directly from filtered.data, one column per thread. */
    ST_double sum_w_c = 0.0;  /* the total of the retained weights, in row order */
    {
        /* rows: the loaded rows retained, or NULL when that is every one of
         * them (no row missing, none a singleton): the columns are then
         * copied as they are */
        const int all_rows = N == N_total;
        ST_int *rows = all_rows ? NULL
                                : (ST_int *)ctools_safe_malloc2((size_t)(N > 0 ? N : 1), sizeof(ST_int));
        if (!all_rows && !rows) {
            SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
        /* fe_sing.keep: N_valid flags, NULL = all kept */
        if (rows && (iv_rc = civ_kept_rows(rows, valid_mask, fe_sing.keep, N_total, N))) {
            SF_error("civreghdfe: could not list the retained rows\n");
            free(rows);
            goto iv_cleanup;
        }
        /* A wrapper that sets __civreghdfe_sample_clear takes its estimation
         * sample from its sample marker (the last plugin variable, 1 on every
         * loaded row): the rows dropped here are cleared in it, one store per
         * dropped row. __civreghdfe_marked tells the wrapper that the marker
         * was updated. */
        {
            ST_double clear_value = 0;
            if (input.touse && SF_scal_use("__civreghdfe_sample_clear", &clear_value) == 0 &&
                clear_value == 1.0) {
                const ST_int marker = SF_nvars();
                ST_int next = 0;
                for (ST_int ii = 0; ii < N_total && !iv_rc && N < N_total; ii++) {
                    if (next < N && rows[next] == ii) { next++; continue; }
                    iv_rc = (_stata_)->safestore(marker, (ST_int)filtered.obs_map[ii], 0.0);
                }
                if (!iv_rc) iv_rc = ctools_scal_save("__civreghdfe_marked", 1.0);
                if (iv_rc) {
                    free(rows);
                    goto iv_cleanup;
                }
            }
        }
        if (rows) {
            /* obs_map in place and in order (a row moves down); the
               clusters into arrays of their own */
            for (ST_int j = 0; j < N; j++) filtered.obs_map[j] = filtered.obs_map[rows[j]];
            if (has_cluster || has_cluster2) {
                ST_int j;
                #pragma omp parallel for schedule(static) if (N > 200000)
                for (j = 0; j < N; j++) {
                    if (has_cluster) cluster_ids_c[j] = cluster_ids[rows[j]];
                    if (has_cluster2) cluster2_ids_c[j] = cluster2_ids[rows[j]];
                }
            }
        } else {
            if (has_cluster) memcpy(cluster_ids_c, cluster_ids, (size_t)N * sizeof(ST_int));
            if (has_cluster2) memcpy(cluster2_ids_c, cluster2_ids, (size_t)N * sizeof(ST_int));
        }
        /* y, X_endog, X_exog, excluded instruments, tsset columns, weights
           (whose total, in row order, the same task adds up) */
        const ST_int n_gather = 1 + K_endog + K_exog + K_excl_in + ts_cols + (has_weights ? 1 : 0);
        ST_int c;
        /* The regression's columns go to their places in all_data, with the
           sums that do not wait for the normalized weights taken on the way
           (civ_column_gather) */
        #pragma omp parallel for schedule(dynamic) if ((size_t)N * (size_t)n_gather > 200000)
        for (c = 0; c < n_gather; c++) {
            ST_int col;
            if (has_weights && c == n_gather - 1) {
                civ_gather(weights_c, filtered.data.vars[var_weight_idx].data.dbl, rows, N);
                ST_double total = 0.0;
                for (ST_int i = 0; i < N; i++) total += weights_c[i];
                sum_w_c = total;
                continue;
            } else if (c == 0) {
                col = col_y;
            } else if (c < 1 + K_endog) {
                col = col_endog + (c - 1);
            } else if (c < 1 + K_endog + K_exog) {
                col = col_exog + (c - 1 - K_endog);
            } else if (c < 1 + K_endog + K_exog + K_excl_in) {
                col = col_excl + (c - 1 - K_endog - K_exog);
            } else {
                ST_int k = c - 1 - K_endog - K_exog - K_excl_in;
                civ_gather(ts_c + (size_t)k * N, filtered.data.vars[var_ts_start_idx + k].data.dbl,
                           rows, N);
                continue;
            }
            /* y, X_endog, X_exog and the excluded instruments are loaded in
               that order from var_y_idx on */
            civ_column *pc = &part_cols[col - col_excl];
            pc->src = filtered.data.vars[var_y_idx + c].data.dbl;
            pc->rows = rows;
            pc->x = all_data + (size_t)col * N;
            /* Z's copy of an exogenous regressor */
            pc->copy = col >= col_exog ? all_data + (size_t)(col - col_exog) * N : NULL;
            pc->N = N;
            civ_column_gather(pc, !has_weights);
            pc->src = NULL;
            pc->rows = NULL;
        }
        if (keep_raw) {
            /* the columns as they are now, for the readers named above */
            const ST_int n_raw = 1 + K_endog + 2 * K_exog + K_excl_in;
            #pragma omp parallel for schedule(dynamic) if ((size_t)N * (size_t)n_raw > 200000)
            for (c = 0; c < n_raw; c++) {
                const ST_double *from;
                ST_double *to;
                if (c == 0) {
                    from = all_data + (size_t)col_y * N;
                    to = y_c;
                } else if (c < 1 + K_endog) {
                    from = all_data + (size_t)(col_endog + c - 1) * N;
                    to = X_endog_c + (size_t)(c - 1) * N;
                } else if (c < 1 + K_endog + K_exog) {
                    from = all_data + (size_t)(col_exog + c - 1 - K_endog) * N;
                    to = X_exog_c + (size_t)(c - 1 - K_endog) * N;
                } else if (c < 1 + K_endog + 2 * K_exog) {
                    /* Z = [X_exog, excluded instruments] */
                    const ST_int k = c - 1 - K_endog - K_exog;
                    from = all_data + (size_t)(col_exog + k) * N;
                    to = Z_c + (size_t)k * N;
                } else {
                    const ST_int k = c - 1 - K_endog - 2 * K_exog;
                    from = all_data + (size_t)(col_excl + k) * N;
                    to = Z_c + (size_t)(K_exog + k) * N;
                }
                memcpy(to, from, (size_t)N * sizeof(ST_double));
            }
        }
        /* slope variables, one array per slope factor */
        for (ST_int g = 0, k = 0; g < G && !iv_rc; g++) {
            if (!input.slope_flag[g]) continue;
            slope_c[g] = (ST_double *)ctools_safe_malloc2((size_t)(N > 0 ? N : 1), sizeof(ST_double));
            if (!slope_c[g]) {
                SF_error("civreghdfe: Memory allocation failed for the slope variables\n");
                iv_rc = 920;
                break;
            }
            civ_gather(slope_c[g], filtered.data.vars[var_slope_start_idx + k++].data.dbl, rows, N);
        }
        free(rows);
        if (iv_rc) goto iv_cleanup;
    }

    /* The FE levels keep their valid rows until the factors are built, which
     * shifts them over the rows dropped here */

    /* Weight normalization:
       - aw (type=1) and pw (type=3): normalize weights so sum(w) = N
       - fw (type=2): compute N_eff = sum(fw) for reporting/DOF */
    ST_double N_eff = (ST_double)N;
    if (has_weights && weights_c) {
        /* sum_w_c: their total in row order, added up during the compaction */
        if (weight_type == 1 || weight_type == 3) {
            /* aweight or pweight: normalize so sum(w) = N */
            const ST_double sum_w_raw = sum_w_c;
            if (sum_w_raw > 0.0) {
                const ST_double scale = (ST_double)N / sum_w_raw;
                ST_int i;
                #pragma omp parallel for schedule(static) if (N > 200000)
                for (i = 0; i < N; i++) weights_c[i] *= scale;
            }
        } else if (weight_type == 2) {
            /* fweight: N_eff = sum(fw) */
            N_eff = sum_w_c;
        }
    }

    /* Retain the compacted observation mapping for residual output. */
    est_obs = filtered.obs_map;
    filtered.obs_map = NULL;

    /* Free original arrays and filtered data (kept alive for direct reads during compaction) */
    free(cluster_ids); free(cluster2_ids); free(valid_mask);
    cluster_ids = cluster2_ids = valid_mask = NULL;
    /* The loaded variables, each as a task of the helper threads while the
       factors are built: handing 560 MB back takes 14 ms in turn here. What
       is left of filtered goes when the tasks are done (before the
       projection, or at the exit). Without the helpers, now. */
    if (filtered.data.vars && filtered.data.nvars > 1 &&
        filtered.data.nobs >= CTOOLS_FREE_THREADED_ROWS && (free_pool = ctools_get_global_pool()) &&
        ctools_persistent_pool_submit_batch(free_pool, civ_variable_free_run, filtered.data.vars,
                                            filtered.data.nvars, sizeof(filtered.data.vars[0])) == 0) {
        /* free_pool: waited for below */
    } else {
        free_pool = NULL;
        ctools_filtered_data_free(&filtered);
    }

    /* Check we still have enough observations */
    if (N <= K_total || N < K_iv) {
        SF_error("civreghdfe: Insufficient observations after singleton removal\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* End singleton timing, start HDFE setup timing */
    t_singleton = ctools_timer_seconds() - t_singleton_start;
    double t_hdfe_setup_start = ctools_timer_seconds();

    /* Initialize the HDFE state on the estimation sample */
    state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (!state) {
        SF_error("civreghdfe: Memory allocation failed for HDFE state\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    state->K = 1 + K_endog + K_exog + K_iv;  /* Total columns to demean */
    state->in1 = 1;
    state->in2 = N;
    state->has_weights = has_weights;
    state->weight_type = weight_type;
    state->weights = weights_c;
    state->maxiter = maxiter;
    state->tolerance = tolerance;
    state->verbose = verbose;

    /* Factors: the CG solver expects levels 1..num_levels, so the levels
       left after singleton removal are renumbered and counted (level_max
       holds each factor's level count before the removal) */
    /* The levels are in order of first appearance among the valid rows: with
       no singleton they are the factors' as they are, with a few they are
       shifted in place. The rows per level are counted here, and the
       weights of the estimation rows (already normalized) summed per level
       in row order */
    HDFE_LoadedFactor fe_loaded[10];
    memset(fe_loaded, 0, sizeof(fe_loaded));
    for (ST_int g = 0; g < G; g++) {
        fe_loaded[g].levels = fe_levels_c[g];
        fe_loaded[g].num_levels = level_max[g];
        fe_loaded[g].num_obs = N_valid;
    }
    ST_int counts_pending = 0;
    ST_retcode factor_rc = hdfe_build_factors_ordered(state, G, fe_loaded, fe_order, &fe_sing, N_valid,
                                                      NULL, has_weights, &counts_pending);
    /* a level array that became a factor's is no longer this function's */
    for (ST_int g = 0; g < G; g++) fe_levels_c[g] = fe_loaded[g].levels;
    if (!factor_rc) {
        if (counts_pending) {
            hdfe_count_factors(state, fe_order, N, NULL, has_weights ? weights_c : NULL);
        } else if (has_weights) {
            for (ST_int g = 0; g < G; g++) hdfe_count_weights(state, g, N, NULL, weights_c);
        }
    }
    free(level_max);
    level_max = NULL;
    for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
    free(fe_levels_c);
    fe_levels_c = NULL;
    hdfe_factor_order_free(fe_order, 10);
    hdfe_singletons_free(&fe_sing);
    if (factor_rc) {
        SF_error("civreghdfe: Memory allocation failed for FE factors\n");
        iv_rc = factor_rc;
        goto iv_cleanup;
    }
    /* Slope factors (before the solver buffers: their counts become the
       sums of squares of the centered slopes per level) */
    if (num_slopes > 0) {
        ST_int *ones = (ST_int *)ctools_safe_malloc2((size_t)(N > 0 ? N : 1), sizeof(ST_int));
        if (!ones) {
            iv_rc = 920;
            goto iv_cleanup;
        }
        for (ST_int i = 0; i < N; i++) ones[i] = 1;
        creghdfe_input ci;
        civ_slope_input(&input, &ci);
        iv_rc = creghdfe_setup_slopes(state, &ci, (const ST_double *const *)slope_c, ones, N,
                                      has_weights ? weights_c : NULL);
        free(ones);
        if (iv_rc) {
            SF_error(iv_rc == 498 ? "civreghdfe: the weights sum to zero\n"
                                  : "civreghdfe: Memory allocation failed for the slope factors\n");
            goto iv_cleanup;
        }
    }

    /* Allocate inv_counts, inv_weighted_counts, and CG solver buffers */
    ST_int num_threads = ctools_get_max_threads();
    if (num_threads < 1) num_threads = 1;
    state->num_threads = num_threads;

    /* Columns partialled out: y, endogenous, exogenous and excluded
       instruments (the buffer count is capped at this number) */
    ST_int max_partial_cols = 1 + K_endog + K_iv;
    if (ctools_hdfe_alloc_buffers(state, 0, max_partial_cols) != 0) {
        SF_error("civreghdfe: Memory allocation failed for HDFE buffers\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    state->factors_initialized = 1;

    /* Record this estimator's active state for lifecycle tracking */
    iv_state = state;

    /* End HDFE setup timing, start FWL timing */
    t_hdfe_setup = ctools_timer_seconds() - t_hdfe_setup_start;
    double t_fwl_start = ctools_timer_seconds();

    /* FWL partialling: partial() arrives as one 0/1 flag per X_exog column
       (matrix __civreghdfe_partial_flags, aligned with the expanded exogenous
       regressors). is_partial: 1 = partialled out, 2 = dropped as redundant. */
    ST_double dval_partial;
    ST_int n_partial = 0;

    if (SF_scal_use("__civreghdfe_n_partial", &dval_partial) == 0) {
        n_partial = (ST_int)dval_partial;
    }

    if (n_partial > 0 && K_exog > 0) {
        is_partial = (ST_int *)calloc(K_exog, sizeof(ST_int));
        ST_int n_flagged = is_partial ? read_flag_matrix("__civreghdfe_partial_flags", K_exog, is_partial) : -1;
        if (n_flagged != n_partial) {
            ST_retcode flag_rc = is_partial ? 198 : 920;
            SF_error(is_partial ? "civreghdfe: partial() column flags do not match the exogenous regressors\n"
                                : "civreghdfe: Memory allocation failed for partial indices\n");
            iv_rc = flag_rc;
            goto iv_cleanup;
        }
        if (bs_retain) {
            if (!(bs_partial = (ST_int *)malloc((size_t)K_exog * sizeof(ST_int)))) {
                iv_rc = 920;
                goto iv_cleanup;
            }
            memcpy(bs_partial, is_partial, (size_t)K_exog * sizeof(ST_int));
        }
    }

    /* End FWL timing, start partial out timing */
    t_fwl = ctools_timer_seconds() - t_fwl_start;
    double t_partial_start = ctools_timer_seconds();

    /* The loaded variables have gone back by now, or do so here */
    if (free_pool) {
        ctools_persistent_pool_wait(free_pool);
        free_pool = NULL;
        ctools_filtered_data_free(&filtered);
    }

    /* Partial out FEs from all variables. The columns are in all_data since
       the compaction; they are weighted as the regression is. */
    for (ST_int p = 0; p < n_part; p++) part_cols[p].w = weights_c;

    /* Pre-partialling (weighted) sums of squares of every column: about the
       mean (ss_pre, for partial()'s redundancy rule, absorb_tol being
       reghdfe's collinear-with-FE tolerance) and uncentred (ss_unc, the
       scale of the rounding noise that the projection leaves in a column
       the fixed effects absorb; see CIV_FE_NOISE) */
    const ST_double absorb_tol = (tolerance / 10.0 < 1e-6) ? tolerance / 10.0 : 1e-6;
    ST_double *ss_pre = all_data + (size_t)N * total_cols;
    ST_double *ss_unc = all_data + (size_t)(N + 1) * total_cols;

    /* Partial out the FEs from every column in one pass, then the partial()
       variables by weighted FWL on the FE-demeaned data: the residual of the
       W-weighted projection on [FE, P] is M_P~ M_FE x with P~ = M_FE P, the
       way ivreghdfe runs HDFE partial_out and then s_partial. */
    /* Each column is projected at unit scale, as chdfe does: a power-of-two
       rescaling keeps the mantissas (the projection is unchanged), but the
       solver's convergence test then behaves the same for a column in units
       of 1e5 as in units of 1. Without it, a column of large values with
       little within variation (1000*id + noise) can come back as noise
       larger than the column itself. */
    HDFE_SolveResult projection = {0, 1, 0};
    ST_double *const part_data = all_data + (size_t)N * col_excl;
    if (state->G == 1 && !state->factors[0].slope && !input.map && state->N == N &&
        state->num_threads >= 1 && maxiter >= 1 && isfinite(tolerance) && tolerance > 0 &&
        !ctools_cancel_requested()) {
        /* A single absorbed effect, which partial_out_columns() projects on
           in one step: each column's remaining passes with that projection
           and the return to its units done in them (civ_column_absorb). The
           level sums are a column's scratch; the whole team runs the loop, a
           column per thread. */
        const FE_Factor *f = &state->factors[0];
        const int team = omp_get_max_threads();
        const size_t lanes = (size_t)(n_part < team ? n_part : team);
        ST_double *level_sums = (ST_double *)ctools_safe_malloc3(lanes, (size_t)f->num_levels,
                                                                  sizeof(ST_double));
        if (!level_sums) {
            SF_error("civreghdfe: Memory allocation failed for demeaning buffer\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
        ST_int p;
        #pragma omp parallel for schedule(static) if (n_part > 1 && (size_t)N * (size_t)n_part > 200000)
        for (p = 0; p < n_part; p++) {
            civ_column *pc = &part_cols[p];
            if (pc->w) civ_column_sums(pc);
            civ_column_absorb(pc, f, level_sums + (size_t)omp_get_thread_num() * (size_t)f->num_levels);
        }
        free(level_sums);
        projection = ctools_cancel_requested() ? (HDFE_SolveResult){1, 0, 1}
                                               : (HDFE_SolveResult){0, 1, 1};
    } else {
        /* A column's passes are one thread's work, so the columns keep n_part
           threads busy. When the projection will be solved for the effects'
           coefficients, what that solve can do before it sees a column
           (hdfe_solver.h) is done meanwhile: the columns run on the helper
           threads and this thread prepares the solve with the OpenMP team. */
        const int ahead = state->G >= 2 && !input.map &&
                          (size_t)N * (size_t)n_part >= CIV_AHEAD_MIN_VALUES;
        ctools_persistent_pool *column_pool = ahead ? ctools_get_global_pool() : NULL;
        if (column_pool &&
            ctools_persistent_pool_submit_batch(column_pool, civ_column_center_run, part_cols,
                                                (size_t)n_part, sizeof(part_cols[0])) == 0) {
            prepared = partial_out_prepare(state, n_part);
            ctools_persistent_pool_wait(column_pool);
        } else {
            ST_int p;
            #pragma omp parallel for schedule(static) if (n_part > 1 && (size_t)N * (size_t)n_part > 200000)
            for (p = 0; p < n_part; p++) civ_column_center_run(&part_cols[p]);
        }
        /* a solver option of reghdfe's given: its default iteration, on the
           rows; otherwise with what was prepared, if anything */
        projection = input.map ? partial_out_columns_options(state, part_data, N, n_part, 0, 1)
                               : partial_out_columns_prepared(state, part_data, N, n_part,
                                                              num_threads, prepared);
        if (prepared) {
            /* its memory goes back after the command */
            ctools_release_later_call(civ_release_prepared, prepared);
            prepared = NULL;
        }
        /* back to the columns' units, and Z's copy of X_exog */
        civ_columns_restore(part_cols, n_part, N);
    }
    for (ST_int p = 0; p < n_part; p++) {
        ss_pre[col_excl + p] = part_cols[p].ss_pre;
        ss_unc[col_excl + p] = part_cols[p].ss_unc;
    }
    for (ST_int k = 0; k < K_exog; k++) {
        ss_pre[k] = ss_pre[col_exog + k];
        ss_unc[k] = ss_unc[col_exog + k];
    }
    const ST_double *ss_unc_endog = ss_unc + col_endog;
    ST_double *ss_exog = ss_pre + col_exog;  /* compacted with X_exog below */
    ST_double *ss_unc_exog = ss_unc + col_exog;  /* likewise */
    const ST_double *ss_unc_excl = ss_unc + col_excl;
    if (bs_retain && !projection.status) {
        const size_t n_tilde = (size_t)N * (size_t)(total_cols - col_excl);
        if (!(bs_tilde = (ST_double *)malloc(n_tilde * sizeof(ST_double)))) {
            iv_rc = 920;
            goto iv_cleanup;
        }
        memcpy(bs_tilde, all_data + (size_t)N * col_excl, n_tilde * sizeof(ST_double));
    }
    /* bw(auto) with absorb(): ivreghdfe selects the bandwidth from the
       instruments after the fixed effects are partialled out (and before
       partial()) */
    if (!projection.status && G > 0 && kernel_type > 0 && bw < 0) {
        abw_fe = (ST_double *)ctools_safe_calloc2((size_t)N, sizeof(ST_double));
        if (!abw_fe) {
            iv_rc = 920;
            goto iv_cleanup;
        }
        for (ST_int c = 0; c < K_excl_in + K_exog; c++) {
            if (c >= K_excl_in && is_partial && is_partial[c - K_excl_in]) continue;
            const ST_double *x = all_data + (size_t)N * (c < K_excl_in ? col_excl + c : col_exog + (c - K_excl_in));
            for (ST_int i = 0; i < N; i++) abw_fe[i] += x[i];
        }
    }
    ST_retcode fwl_rc = 0;
    if (!projection.status && is_partial) {
        fwl_rc = civ_fwl_partial_out(all_data, N, total_cols, col_exog, 0, K_exog,
                                     is_partial, weights_c, ss_exog, absorb_tol, 1);
    }

    if (projection.status || fwl_rc) {
        if (projection.status && projection.status != 1) SF_error("civreghdfe: fixed-effect projection did not converge\n");
        iv_rc = projection.status ? projection.status : fwl_rc;
        goto iv_cleanup;
    }

    /* End partial out timing, start post-processing timing */
    t_partial_out = ctools_timer_seconds() - t_partial_start;
    double t_postproc_start = ctools_timer_seconds();

    /* Extract demeaned data */
    ST_double *y_dem = all_data + (size_t)N * col_y;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + (size_t)N * col_endog : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * col_exog : NULL;
    ST_double *Z_dem = all_data;

    /* Remove partial variables from X_exog AND from Z (instruments) after convergence.
       Z contains [exog_vars, excluded_instruments], so partial exog vars must be removed
       from both X_exog and the first K_exog columns of Z. Partial columns dropped as
       redundant (flag 2) are removed as well, and reported back for the ado's note. */
    if (is_partial) {
        ST_int K_exog_new = K_exog - n_partial;
        ST_int K_excl = K_iv - K_exog;  /* Number of excluded instruments */

        for (ST_int k = 0; k < K_exog; k++) {
            if (is_partial[k] == 2 &&
                (iv_rc = ctools_mat_store("__civreghdfe_partial_flags", 1, k + 1, 2.0)))
                goto iv_cleanup;
        }

        /* Remove partial columns from X_exog (and their pre-partialling SS) */
        if (K_exog_new > 0) {
            ST_int new_idx = 0;
            for (ST_int k = 0; k < K_exog; k++) {
                if (!is_partial[k]) {
                    if (new_idx != k) {
                        memcpy(X_exog_dem + (size_t)new_idx * N, X_exog_dem + (size_t)k * N, N * sizeof(ST_double));
                        ss_exog[new_idx] = ss_exog[k];
                        ss_unc_exog[new_idx] = ss_unc_exog[k];
                    }
                    new_idx++;
                }
            }
        }

        /* Remove partial columns from Z (which starts with K_exog columns of exog vars) */
        /* After removal, Z should have K_exog_new + K_excl columns */
        ST_int new_idx = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_partial[k]) {
                if (new_idx != k) {
                    memcpy(Z_dem + (size_t)new_idx * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_idx++;
            }
        }
        /* Shift excluded instruments down */
        for (ST_int k = 0; k < K_excl; k++) {
            memcpy(Z_dem + (size_t)new_idx * N, Z_dem + (size_t)(K_exog + k) * N, N * sizeof(ST_double));
            new_idx++;
        }

        K_exog = K_exog_new;
        K_iv = K_exog_new + K_excl;
        X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * col_exog : NULL;
        /* is_partial stays until the absorbed columns are counted below: it
           maps the remaining X_exog columns back to X_exog_c */
    }

    /* Compute df_a (absorbed degrees of freedom) and mobility groups. The
       redundant coefficients of each effect are counted as reghdfe counts
       them (the HDFE object behind ivreghdfe's absorb_ct), with every pair of
       effects examined; the counts are redone below once the cluster nesting
       is known. */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    ST_int mobility_groups = 0;

    /* dofadjustments(), as reghdfe (which ivreghdfe calls) applies them:
       0 = all (pairwise and clusters), 1 = none, 2 = firstpair,
       3 = pairwise; dof_clusters: effects nested in the clusters */
    ST_double dof_type_d = 0, dof_clusters_d = 1;
    SF_scal_use("__civreghdfe_dof_type", &dof_type_d);
    SF_scal_use("__civreghdfe_dof_clusters", &dof_clusters_d);
    const ST_int dof_pairs = dof_type_d == 1 ? 0 :
        (dof_type_d == 2 ? CTOOLS_DOF_FIRSTPAIR : CTOOLS_DOF_PAIRWISE);
    const ST_int dof_clusters = dof_clusters_d != 0;

    fe_redundant = (ST_int *)ctools_safe_calloc2((size_t)(G > 0 ? G : 1), sizeof(ST_int));
    fe_nested = (ST_int *)ctools_safe_calloc2((size_t)(G > 0 ? G : 1), sizeof(ST_int));
    fe_exact = (ST_int *)ctools_safe_calloc2((size_t)(G > 0 ? G : 1), sizeof(ST_int));
    if (!fe_redundant || !fe_nested || !fe_exact ||
        civ_redundant(state, &input, NULL, dof_pairs, fe_redundant, fe_exact)) {
        SF_error("civreghdfe: Memory allocation failed for the absorbed degrees of freedom\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    for (ST_int g = 0; g < G; g++) df_a += state->factors[g].num_levels - fe_redundant[g];
    {
        /* the mobility groups of the second intercept factor, and one
           redundant coefficient for each further one (pairwise) */
        ST_int n_int = 0, second = -1;
        for (ST_int g = 0; g < G; g++) {
            if (input.slope_flag[g]) continue;
            if (++n_int == 2) second = g;
        }
        if (second >= 0 && dof_pairs)
            mobility_groups = fe_redundant[second] + (dof_pairs == CTOOLS_DOF_PAIRWISE ? n_int - 2 : 0);
    }

    for (ST_int g = 0; g < G; g++) {
        /* Save per-FE num_levels to Stata scalars for absorbed DOF table */
        char scalar_name[64];
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_num_levels_%d", (int)(g + 1));
        if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)state->factors[g].num_levels)))
            goto iv_cleanup;
    }

    /* Save mobility groups for absorbed DOF table display */
    if ((iv_rc = ctools_scal_save("__civreghdfe_mobility_groups", (ST_double)mobility_groups)))
        goto iv_cleanup;

    /* For VCE calculation when FE is nested in cluster:
       - nested FE levels don't contribute to df_a for VCE
       - nested_adj = 1 to account for the constant */
    ST_int df_a_for_vce = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;

    /* Remap cluster IDs to contiguous 1-based indices using shared utility */
    ST_int num_clusters = 0;
    if (has_cluster) {
        if (ctools_remap_cluster_ids(cluster_ids_c, N, &num_clusters) != 0) {
            SF_error("civreghdfe: Failed to remap cluster IDs (memory allocation error)\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
    }

    /* Remap cluster2 IDs for two-way clustering */
    ST_int num_clusters2 = 0;
    if (has_cluster2) {
        if (ctools_remap_cluster_ids(cluster2_ids_c, N, &num_clusters2) != 0) {
            SF_error("civreghdfe: Failed to remap cluster2 IDs (memory allocation error)\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
    }

    if (!kiefer && (vce_type == 2 || vce_type == 3) &&
        (num_clusters < 2 || (has_cluster2 && num_clusters2 < 2))) {
        SF_error("civreghdfe: clustered VCE requires at least two retained clusters in each dimension\n");
        iv_rc = 459;
        goto iv_cleanup;
    }

    /* Detect FEs nested within cluster variable using data-based check.
     * This works for both numeric and string cluster variables. As in
     * reghdfe, an effect nested in either dimension of a two-way cluster is
     * wholly redundant, and the redundant coefficients of the other effects
     * are counted with the nested ones set aside. df_a keeps the nested
     * levels; df_a_for_vce = df_a - df_a_nested is reghdfe's df_a. Kiefer's
     * "cluster" is the panel variable of the HAC lag structure, not a
     * clustering, so nothing is nested there (ivreghdfe's absorb_ct is the
     * full df_a). */
    /* A failed allocation must not silently skip the detection: it changes
     * the cluster dof */
    if (has_cluster && !kiefer && dof_clusters) {
        df_a_nested = 0;
        for (ST_int g = 0; g < G; g++) {
            /* a slope is never nested in the clusters */
            if (input.slope_flag[g]) {
                fe_nested[g] = 0;
                continue;
            }
            ST_int is_nested = ctools_fe_nested_in_cluster(
                state->factors[g].levels, state->factors[g].num_levels,
                cluster_ids_c, N);
            if (is_nested == 0 && has_cluster2)
                is_nested = ctools_fe_nested_in_cluster(
                    state->factors[g].levels, state->factors[g].num_levels,
                    cluster2_ids_c, N);
            if (is_nested < 0) {
                SF_error("civreghdfe: Memory allocation failed for the nesting check\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
            fe_nested[g] = is_nested;
            if (is_nested) df_a_nested += state->factors[g].num_levels;
        }
        if (df_a_nested > 0) {
            if (civ_redundant(state, &input, fe_nested, dof_pairs, fe_redundant, fe_exact)) {
                SF_error("civreghdfe: Memory allocation failed for the absorbed degrees of freedom\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
            df_a = 0;
            for (ST_int g = 0; g < G; g++)
                df_a += state->factors[g].num_levels - (fe_nested[g] ? 0 : fe_redundant[g]);
        }
        /* Recompute VCE DOF with data-based nesting results */
        df_a_for_vce = df_a - df_a_nested;
        nested_adj = (df_a_nested > 0) ? 1 : 0;
    }
    for (ST_int g = 0; g < G; g++) {
        char scalar_name[64];
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_nested_%d", (int)(g + 1));
        if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)fe_nested[g])))
            goto iv_cleanup;
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_redundant_%d", (int)(g + 1));
        if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)fe_redundant[g])))
            goto iv_cleanup;
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_exact_%d", (int)(g + 1));
        if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)fe_exact[g])))
            goto iv_cleanup;
    }

    /* For Kiefer, the panel variable is set as "cluster" for HAC lag structure,
       not true clustering. The FE nesting detection incorrectly sets df_a_for_vce=0
       because the FE and "cluster" are the same variable. Override: use full df_a
       to match ivreghdfe's sdofminus = absorb_ct. */
    if (kiefer) {
        df_a_for_vce = df_a;
        nested_adj = 0;
    }

    /* Update K_total after possible FWL reduction of K_exog */
    K_total = K_exog + K_endog;

    /* ================================================================
     * COLLINEARITY DETECTION for X (regressors) and Z (instruments)
     * Detect after HDFE partialling, before estimation.
     * Same approach as creghdfe: FE-absorbed variance check + Cholesky.
     * ================================================================ */
    ST_int K_exog_orig = K_exog;
    ST_int K_endog_orig = K_endog;
    ST_int K_total_orig = K_total;

    /* is_collinear_x: flags for [exog, endog] in internal order */
    is_collinear_x = (ST_int *)calloc(K_total, sizeof(ST_int));
    ST_int num_collinear_x = 0;

    if (!is_collinear_x) {
        SF_error("civreghdfe: Memory allocation failed for collinearity arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Stage 1: columns absorbed by the fixed effects (or by partial()): the
       partialled (weighted) SS is rounding noise relative to the column's
       uncentred SS (CIV_FE_NOISE). They are omitted and do not count in the
       ranks, as ivreghdfe's invsym() zeroes them (when its own projection
       leaves noise above invsym's tolerance, ivreghdfe keeps a column of
       noise with an arbitrary coefficient; the model without the column is
       the reference). Every absorbed or collinear exogenous regressor still
       counts in the Cragg-Donald and KP Wald F dof, whose instrument count
       ivreghdfe takes from its variable lists (ranktest partials those
       regressors out, so the statistics themselves do not change). An
       absorbed excluded instrument makes ivreghdfe's statistics missing
       (exact zero), so there the model without it is the reference. */
    ST_int n_abs_x_exog = 0, n_abs_x_endog = 0, n_abs_z_excl = 0, n_abs_exog_all = 0;
    for (ST_int k = 0; k < K_exog; k++) {
        ST_double xx_partial = civ_wss(X_exog_dem + (size_t)k * N, weights_c, N);
        if (!(xx_partial > CIV_FE_NOISE * ss_unc_exog[k])) {
            is_collinear_x[k] = 1;
            num_collinear_x++;
        }
    }
    for (ST_int k = 0; k < K_endog; k++) {
        ST_double xx_partial = civ_wss(X_endog_dem + (size_t)k * N, weights_c, N);
        if (!(xx_partial > CIV_FE_NOISE * ss_unc_endog[k])) {
            is_collinear_x[K_exog + k] = 1;
            num_collinear_x++;
        }
    }
    free(is_partial);
    is_partial = NULL;

    /* Stage 2: ivreghdfe's rank rule, invsym() of the partialled X'WX in its
       column order [endogenous, exogenous] (the endogenous block precedes
       the exogenous one in all_data), with the stage-1 columns left out.
       Exact collinearity among the variables themselves is resolved by the
       ado first, in ivreg2's order; what remains here is collinearity that
       the fixed effects create. */
    if (K_total > 0) {
        ST_double *XtX = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
        ST_int *drop = (ST_int *)calloc((size_t)K_total, sizeof(ST_int));
        ST_int n_drop = -1;
        if (XtX && drop) {
            const ST_double *X_ne = (K_endog > 0) ? X_endog_dem : X_exog_dem;
            ctools_matmul_atdb(X_ne, X_ne, weights_c, N, K_total, K_total, XtX);
            for (ST_int k = 0; k < K_endog; k++) drop[k] = is_collinear_x[K_exog + k];
            for (ST_int k = 0; k < K_exog; k++) drop[K_endog + k] = is_collinear_x[k];
            n_drop = civ_invsym_drops(XtX, K_total, 0, drop);
            for (ST_int k = 0; k < K_endog; k++) is_collinear_x[K_exog + k] = drop[k];
            for (ST_int k = 0; k < K_exog; k++) is_collinear_x[k] = drop[K_endog + k];
        }
        free(XtX);
        free(drop);
        if (n_drop < 0) {
            SF_error("civreghdfe: Memory allocation failed for the collinearity check\n");
            iv_rc = 920;
            goto iv_cleanup;
        }

        /* Recount total collinear */
        num_collinear_x = 0;
        for (ST_int k = 0; k < K_total; k++) {
            if (is_collinear_x[k]) num_collinear_x++;
        }
        for (ST_int k = 0; k < K_exog; k++) n_abs_exog_all += is_collinear_x[k];
    }

    /* Store collinearity flags to Stata scalars (in [exog, endog] order) */
    {
        char scalar_name[64];
        if ((iv_rc = ctools_scal_save("__civreghdfe_num_collinear", (ST_double)num_collinear_x)))
            goto iv_cleanup;
        for (ST_int k = 0; k < K_total_orig; k++) {
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_collinear_%d", (int)(k + 1));
            if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)is_collinear_x[k])))
                goto iv_cleanup;
        }
    }

    /* Compact X arrays in-place if collinear columns found */
    if (num_collinear_x > 0) {
        /* Compact X_exog_dem */
        ST_int new_K_exog = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_collinear_x[k]) {
                if (new_K_exog != k) {
                    memcpy(X_exog_dem + (size_t)new_K_exog * N, X_exog_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_K_exog++;
            }
        }

        /* Compact X_endog_dem */
        ST_int new_K_endog = 0;
        for (ST_int k = 0; k < K_endog; k++) {
            if (!is_collinear_x[K_exog + k]) {
                if (new_K_endog != k) {
                    memcpy(X_endog_dem + (size_t)new_K_endog * N, X_endog_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_K_endog++;
            }
        }

        /* Compact exogenous portion of Z_dem (first K_exog columns match X_exog) */
        ST_int K_excl = K_iv - K_exog;
        ST_int new_z_idx = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_collinear_x[k]) {
                if (new_z_idx != k) {
                    memcpy(Z_dem + (size_t)new_z_idx * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_z_idx++;
            }
        }
        /* Shift excluded instruments down after compacted exog portion */
        for (ST_int k = 0; k < K_excl; k++) {
            if (new_z_idx != K_exog + k) {
                memcpy(Z_dem + (size_t)new_z_idx * N, Z_dem + (size_t)(K_exog + k) * N, N * sizeof(ST_double));
            }
            new_z_idx++;
        }

        K_exog = new_K_exog;
        K_endog = new_K_endog;
        K_total = K_exog + K_endog;
        K_iv = new_z_idx;  /* new_K_exog + K_excl */

        if (verbose) {
            char msg[128];
            snprintf(msg, sizeof(msg), "civreghdfe: Dropped %d collinear regressor(s)\n",
                     (int)num_collinear_x);
            SF_display(msg);
        }
    }

    /* Z collinearity detection: excluded instruments absorbed by the fixed
       effects (the stage-1 rule above), then ivreghdfe's invsym() rule on
       Z'WZ. The exogenous block of Z was checked with X. */
    ST_int K_excl_z = K_iv - K_exog;  /* excluded instruments as sent by the ado */
    ST_int K_excl_dropped = 0;
    ST_int *is_collinear_z = (ST_int *)calloc((size_t)(K_iv > 0 ? K_iv : 1), sizeof(ST_int));
    ST_int num_collinear_z = 0;

    if (is_collinear_z) {
        for (ST_int j = 0; j < K_excl_z; j++) {
            ST_double zz_partial = civ_wss(Z_dem + (size_t)(K_exog + j) * N, weights_c, N);
            if (!(zz_partial > CIV_FE_NOISE * ss_unc_excl[j])) is_collinear_z[K_exog + j] = 1;
        }

        /* invsym() of Z'WZ: the kept exogenous regressors are swept first,
           so only excluded instruments are dropped */
        ST_double *ZtZ = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
        ST_int z_drops = ZtZ ? 0 : -1;
        if (ZtZ) {
            ctools_matmul_atdb(Z_dem, Z_dem, weights_c, N, K_iv, K_iv, ZtZ);
            z_drops = civ_invsym_drops(ZtZ, K_iv, K_exog, is_collinear_z);
        }
        free(ZtZ);
        if (z_drops < 0) {
            free(is_collinear_z);
            SF_error("civreghdfe: Memory allocation failed for the collinearity check\n");
            iv_rc = 920;
            goto iv_cleanup;
        }

        /* Compact Z_dem in-place */
        ST_int new_K_iv = 0;
        for (ST_int k = 0; k < K_iv; k++) {
            if (!is_collinear_z[k]) {
                if (new_K_iv != k) {
                    memcpy(Z_dem + (size_t)new_K_iv * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_K_iv++;
            }
        }
        num_collinear_z = K_iv - new_K_iv;
        K_excl_dropped = num_collinear_z;
        /* which excluded instruments were dropped (smatrix(), sfirst names) */
        for (ST_int j = 0; j < K_excl_z && !iv_rc; j++) {
            char zn[64];
            snprintf(zn, sizeof(zn), "__civreghdfe_zcollin_%d", (int)(j + 1));
            iv_rc = ctools_scal_save(zn, (ST_double)(is_collinear_z[K_exog + j] != 0));
        }
        if (iv_rc) {
            free(is_collinear_z);
            goto iv_cleanup;
        }
        K_iv = new_K_iv;
        /* excluded instruments dropped here, for ivreghdfe's footer note */
        if ((iv_rc = ctools_scal_save("__civreghdfe_z_dropped", (ST_double)num_collinear_z))) {
            free(is_collinear_z);
            goto iv_cleanup;
        }

        if (verbose && num_collinear_z > 0) {
            char msg[128];
            snprintf(msg, sizeof(msg), "civreghdfe: Dropped %d collinear instrument(s)\n",
                     (int)num_collinear_z);
            SF_display(msg);
        }
    }

    /* An endogenous regressor collinear with the instruments, which ivreg2
       reclassifies as exogenous when the variables themselves are collinear:
       flagged for the ado, which then applies ivreg2's rules (invsym() swept
       in order over [instruments, endogenous]) */
    if (K_endog > 0 && K_iv > 0) {
        ST_int KJ = K_iv + K_endog;
        ST_double *J = (ST_double *)ctools_safe_calloc3((size_t)KJ, (size_t)KJ, sizeof(ST_double));
        ST_double *blk_zx = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_endog, sizeof(ST_double));
        ST_double *blk_xx = (ST_double *)ctools_safe_calloc3((size_t)K_endog, (size_t)K_endog, sizeof(ST_double));
        ST_int *jdrop = (ST_int *)calloc((size_t)KJ, sizeof(ST_int));
        ST_int j_drops = -1;
        if (J && blk_zx && blk_xx && jdrop) {
            ST_double *blk_zz = J;  /* filled in place below */
            ST_double *zz = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
            if (zz) {
                ctools_matmul_atdb(Z_dem, Z_dem, weights_c, N, K_iv, K_iv, zz);
                ctools_matmul_atdb(Z_dem, X_endog_dem, weights_c, N, K_iv, K_endog, blk_zx);
                ctools_matmul_atdb(X_endog_dem, X_endog_dem, weights_c, N, K_endog, K_endog, blk_xx);
                for (ST_int c = 0; c < KJ; c++)
                    for (ST_int r = 0; r < KJ; r++) {
                        ST_double v;
                        if (c < K_iv && r < K_iv) v = zz[(size_t)c * K_iv + r];
                        else if (c >= K_iv && r >= K_iv) v = blk_xx[(size_t)(c - K_iv) * K_endog + (r - K_iv)];
                        else if (c >= K_iv) v = blk_zx[(size_t)(c - K_iv) * K_iv + r];
                        else v = blk_zx[(size_t)(r - K_iv) * K_iv + c];
                        blk_zz[(size_t)c * KJ + r] = v;
                    }
                j_drops = civ_invsym_drops(J, KJ, KJ, jdrop);
                free(zz);
            }
        }
        free(J); free(blk_zx); free(blk_xx); free(jdrop);
        if (j_drops < 0) {
            free(is_collinear_z);
            SF_error("civreghdfe: Memory allocation failed for the collinearity check\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
        bs_endo_in_z = j_drops;
        if ((iv_rc = ctools_scal_save("__civreghdfe_endo_in_z", (ST_double)j_drops))) {
            free(is_collinear_z);
            goto iv_cleanup;
        }
    }

    /* orthog()/endogtest()/redundant(): map the ado's per-column flags onto
       the excluded instruments and endogenous regressors that survived the
       checks above; test_idx holds the three position lists back to back. */
    ST_int n_orthog_t = 0, n_endogtest_t = 0, n_redund_t = 0;
    {
        ST_double dv;
        if (SF_scal_use("__civreghdfe_n_orthog", &dv) == 0 && dv > 0) n_orthog_t = (ST_int)dv;
        if (SF_scal_use("__civreghdfe_n_endogtest", &dv) == 0 && dv > 0) n_endogtest_t = (ST_int)dv;
        if (SF_scal_use("__civreghdfe_n_redundant", &dv) == 0 && dv > 0) n_redund_t = (ST_int)dv;
    }
    test_idx = (ST_int *)malloc(((size_t)n_orthog_t + n_endogtest_t + n_redund_t + 1) * sizeof(ST_int));
    ST_retcode test_rc = test_idx ? 0 : 920;
    if (test_rc) SF_error("civreghdfe: Memory allocation failed for test columns\n");
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_orthog_flags", "orthog", 1, n_orthog_t, K_excl_z,
                                is_collinear_z ? is_collinear_z + K_exog : NULL, test_idx);
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_endogtest_flags", "endogtest", 2, n_endogtest_t, K_endog_orig,
                                is_collinear_x + K_exog_orig, test_idx + n_orthog_t);
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_redundant_flags", "redundant", 3, n_redund_t, K_excl_z,
                                is_collinear_z ? is_collinear_z + K_exog : NULL,
                                test_idx + n_orthog_t + n_endogtest_t);
    if (bs_retain && !test_rc) {
        bs_collin_z = (ST_int *)calloc((size_t)(K_excl_z > 0 ? K_excl_z : 1), sizeof(ST_int));
        if (!bs_collin_z) test_rc = 920;
        for (ST_int j = 0; bs_collin_z && is_collinear_z && j < K_excl_z; j++)
            bs_collin_z[j] = is_collinear_z[K_exog + j] != 0;
    }
    free(is_collinear_z);

    /* Post-compaction identification check */
    if (K_iv < K_total || test_rc) {
        if (!test_rc) SF_error("civreghdfe: Model is underidentified after removing collinear variables\n");
        iv_rc = test_rc ? test_rc : 481;
        goto iv_cleanup;
    }

    /* Store compacted dimensions; ivreg2's rankzz is the instrument count */
    if ((iv_rc = ctools_scal_save("__civreghdfe_K_keep", (ST_double)K_total)) ||
        (iv_rc = ctools_scal_save("__civreghdfe_rankzz", (ST_double)K_iv)))
        goto iv_cleanup;

    /* Allocate output arrays */
    beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    V = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
    first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Check allocations - critical for preventing NULL pointer dereference */
    if (!beta || !V || !first_stage_F) {
        SF_error("civreghdfe: Memory allocation failed for output arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* End post-processing timing, start estimation timing */
    t_postproc = ctools_timer_seconds() - t_postproc_start;
    double t_estimate_start = ctools_timer_seconds();

    /* Compute k-class IV estimation (2SLS, LIML, Fuller, etc.) */
    /* For cluster VCE, pass df_a_for_vce (excluding nested FE levels) and nested_adj */
    ST_retcode rc = 0;

    /* tsset structure for the kernel estimators: panel ids and the time
       index (t - tmin)/tdelta, so that kernel lags are time differences
       within a panel (ivreg2's L. operators) and Driscoll-Kraay periods
       follow calendar order */
    civreghdfe_ts ts_info;
    memset(&ts_info, 0, sizeof(ts_info));
    ts_int = NULL;  /* panel ids (N), then time index (N) */
    if (ts_mode > 0) {
        ts_int = (ST_int *)ctools_safe_malloc3((size_t)N, 2, sizeof(ST_int));
        if (!ts_int) {
            SF_error("civreghdfe: Memory allocation failed for the time structure\n");
            rc = 920;
        } else {
            const ST_double *tv = ts_c + (size_t)(ts_cols - 1) * N;
            ST_int *tidx = ts_int + N;
            ST_double tmin = tv[0], tmax = tv[0];
            for (ST_int i = 1; i < N; i++) {
                if (tv[i] < tmin) tmin = tv[i];
                if (tv[i] > tmax) tmax = tv[i];
            }
            /* The time index is an int: refuse spans of 2^31-2 or more delta
               units, which would wrap it and break the kernel lag pairing
               (e.g. a %tc clock with delta 1 ms over more than 24.8 days) */
            if (!((tmax - tmin) / tdelta < 2147483646.0)) {
                SF_error("civreghdfe: the tsset time span is too long for kernel estimators "
                         "(2^31-2 or more periods of delta); tsset with a coarser delta()\n");
                rc = 198;
            }
            ST_int tmax_idx = 0;
            for (ST_int i = 0; i < N && !rc; i++) {
                tidx[i] = (ST_int)llround((tv[i] - tmin) / tdelta);
                if (tidx[i] > tmax_idx) tmax_idx = tidx[i];
            }
            /* Time span of the estimation sample: the ado reports T = span + 1
               as Kiefer's bandwidth; as ivreg2, a bandwidth may not exceed
               span / delta */
            if (!rc) rc = ctools_scal_save("__civreghdfe_tspan", tmax - tmin);
            if (!rc && kernel_type > 0 && !kiefer && bw > (tmax - tmin) / tdelta) {
                SF_error("invalid bandwidth in option bw() - cannot exceed timespan of data\n");
                rc = 198;
            }
            ST_int num_panels = 1;
            if (ts_mode == 2 && !rc) {
                if (ctools_numeric_to_cluster_ids(ts_c, N, ts_int, &num_panels) != 0) {
                    SF_error("civreghdfe: Memory allocation failed for the panel structure\n");
                    rc = 920;
                }
                for (ST_int i = 0; i < N && !rc; i++) ts_int[i] += 1;
            }
            ts_info.panel = (ts_mode == 2) ? ts_int : NULL;
            ts_info.num_panels = num_panels;
            ts_info.tidx = tidx;
            ts_info.tmax_idx = tmax_idx;
            ts_info.tdelta = tdelta;
        }
    }

    civreghdfe_test_cols test_cols = {
        test_idx, n_orthog_t,
        test_idx + n_orthog_t, n_endogtest_t,
        test_idx + n_orthog_t + n_endogtest_t, n_redund_t,
        NULL, 0
    };

    /* ivreghdfe with absorb() implies small and partialsmall: its small-sample
       sdofminus is sdofminus() + the partialled-out regressors + absorb_ct
       (max(1, df_a), the nested-in-cluster convention), and nopartialsmall
       is ignored. dofminus() is the large-sample adjustment. */
    ST_int sdof_extra = sdofminus_opt + ((G > 0 || !nopartialsmall) ? n_partial : 0);
    ST_double small_value = 0, has_cons_value = 0;
    SF_scal_use("__civreghdfe_small", &small_value);
    SF_scal_use("__civreghdfe_has_cons", &has_cons_value);
    ST_int small = G > 0 || small_value != 0;
    ST_int dofminus_vce = dofminus;
    ST_int df_a_small = df_a_for_vce + sdof_extra;
    if (sdof_extra > 0 && df_a_for_vce == 0 && nested_adj == 1) df_a_small = 1 + sdof_extra;

    /* b0(): the coefficients at which ivreg2 evaluates the CUE objective,
       in e(b)'s [endogenous, exogenous] order of the plugin columns; the
       estimation works on the kept columns in [exogenous, endogenous] order */
    ST_double *b0 = NULL;
    {
        ST_double has_b0 = 0;
        if (!rc && SF_scal_use("__civreghdfe_has_b0", &has_b0) == 0 && has_b0 == 1) {
            b0 = (ST_double *)calloc((size_t)(K_total > 0 ? K_total : 1), sizeof(ST_double));
            if (!b0) rc = 920;
            ST_int ce = 0, cx = 0;
            for (ST_int i = 0; !rc && i < K_endog_orig + K_exog_orig; i++) {
                ST_int is_endo = i < K_endog_orig;
                ST_int fi = is_endo ? K_exog_orig + i : i - K_endog_orig;
                if (is_collinear_x[fi]) continue;
                ST_double v = 0;
                if (SF_mat_el("__civreghdfe_b0", 1, i + 1, &v)) {
                    SF_error("civreghdfe: b0() could not be read\n");
                    rc = 198;
                    break;
                }
                if (is_endo) b0[K_exog + ce++] = v;
                else b0[cx++] = v;
            }
        }
    }

    civ_extra_buf xb;
    memset(&xb, 0, sizeof(xb));
    if (!rc) rc = civ_build_extra(&xb, N, input.K_exog, K_excl_in, keep_raw ? Z_c : NULL, G, K_exog_orig,
                                  is_collinear_x, K_iv, kernel_type, bw, abw_fe);
    test_cols.orthog_x = xb.ox;
    test_cols.n_orthog_x = xb.n_ox;

    ST_double lambda = 1.0;
    if (!rc) rc = ivest_compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, N_eff, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        cluster2_ids_c, num_clusters2,
        df_a_small, nested_adj, verbose,
        est_method, kclass_user, fuller_alpha, &lambda,
        kernel_type, bw, kiefer,
        ts_mode > 0 ? &ts_info : NULL,
        (G > 0) ? (df_a_for_vce > 0 ? df_a_for_vce : 1) + sdof_extra : sdofminus_opt, center,
        &test_cols,
        n_abs_x_exog + n_abs_x_endog, n_abs_x_exog + n_abs_z_excl, n_abs_exog_all,
        dofminus_vce, small, b0, K_excl_dropped, &xb.ex
    );
    civ_extra_free(&xb);
    free(b0);
    free(test_idx);
    free(ts_int);
    test_idx = ts_int = NULL;

    if (rc != STATA_OK) {
        iv_rc = rc;
        goto iv_cleanup;
    }

    /* An endogenous regressor dropped here (absorbed by the fixed effects,
       or collinear after partialling) has no first-stage variation: the
       rank condition fails, so the identification statistics are 0, as
       ivreghdfe's ranktest finds them with the zero column (up to rounding),
       and their dof count it */
    {
        ST_int n_endo_drop = 0;
        for (ST_int k = 0; k < K_endog_orig; k++) n_endo_drop += is_collinear_x[K_exog_orig + k] != 0;
        if (n_endo_drop > 0) {
            ST_int id_df = (K_iv - K_exog + K_excl_dropped) - K_endog_orig + 1;
            if (id_df < 1) id_df = 1;
            if ((iv_rc = ctools_scal_save("__civreghdfe_underid", 0.0)) ||
                (iv_rc = ctools_scal_save("__civreghdfe_underid_df", (ST_double)id_df)) ||
                (iv_rc = ctools_scal_save("__civreghdfe_cd_f", 0.0)) ||
                (iv_rc = ctools_scal_save("__civreghdfe_kp_f", 0.0)))
                goto iv_cleanup;
        }
    }

    /* End estimation timing (includes VCE), start stats computation */
    t_estimate = ctools_timer_seconds() - t_estimate_start;
    double t_stats_start = ctools_timer_seconds();

    /* Compute RSS, TSS, R-squared, Root MSE, and F-statistic */
    /* These are computed on demeaned (partialled-out) data */
    ST_double rss = 0.0;
    ST_double tss = 0.0;
    ST_double y_mean = 0.0;
    ST_double sum_w = 0.0;

    /* Compute weighted mean of y_dem */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        y_mean += w * y_dem[i];
        sum_w += w;
    }
    y_mean = (G > 0 || has_cons_value != 0) ? y_mean / sum_w : 0.0;

    /* Compute TSS (centered) = sum(w * (y - ybar)^2) */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        ST_double dev = y_dem[i] - y_mean;
        tss += w * dev * dev;
    }

    ST_double resid_idx_value = 0;
    SF_scal_use("__civreghdfe_resid_idx", &resid_idx_value);
    ST_int resid_idx = (ST_int)resid_idx_value;
    ST_retcode residual_rc = 0;
    /* Residuals y - (X_exog * beta[0:K_exog-1] + X_endog * beta[K_exog:K_total-1]),
       each fitted value accumulated in that column order, then
       RSS = sum(w * resid^2) in row order */
    ST_double *fit = (ST_double *)ctools_safe_calloc2((size_t)N, sizeof(ST_double));
    if (!fit) {
        SF_error("civreghdfe: Memory allocation failed for residuals\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    if (K_exog > 0) ctools_matvec_add(X_exog_dem, N, K_exog, beta, fit);
    if (K_endog > 0) ctools_matvec_add(X_endog_dem, N, K_endog, beta + K_exog, fit);
    {
        ST_int row;
        #pragma omp parallel for schedule(static) if (N > 200000)
        for (row = 0; row < N; row++) fit[row] = y_dem[row] - fit[row];
    }
    /* The rows of the estimation sample are distinct observations: stored on
     * all threads (ctools_store_filtered_rowpar proves that first) */
    if (resid_idx > 0)
        residual_rc = ctools_stata_rc(
            ctools_store_filtered_rowpar(fit, (size_t)N, resid_idx, est_obs));
    for (ST_int i = 0; i < N; i++) {
        ST_double resid = fit[i];
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        rss += w * resid * resid;
    }
    free(fit);

    ST_double r2 = (tss > 0) ? 1.0 - rss / tss : 0.0;

    /* Compute df_r:
       - For clustered VCE (vce_type == 2, 3) and Driscoll-Kraay (4, clusters
         = the time periods of the estimation sample): df_r = num_clusters - 1
       - Otherwise: df_r = N - K - df_a - dofminus
       - nopartialsmall: exclude partialled variables from K for DOF calc */
    /* FE-absorbed columns kept as noise by ivreghdfe (stage 1 above) count */
    ST_int n_abs_x = n_abs_x_exog + n_abs_x_endog;
    ST_int K_for_dof = (G > 0) ? K_total + n_abs_x + sdof_extra
                               : K_total + n_abs_x + sdof_extra;
    /* In double: N_eff (the sum of fweights) may exceed INT_MAX */
    ST_double df_r_val;
    if (kiefer) {
        /* Kiefer VCE: df_r = N_eff - K - df_a - dofminus (not cluster-based) */
        df_r_val = N_eff - K_for_dof - df_a - dofminus;
    } else if (has_cluster2) {
        /* Two-way clustering: df_r = min(G1, G2) - 1 */
        ST_int min_clust = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
        df_r_val = min_clust - 1;
    } else if (has_cluster) {
        df_r_val = num_clusters - 1;
    } else {
        df_r_val = N_eff - K_for_dof - df_a - dofminus;
    }
    if (df_r_val <= 0) df_r_val = 1;
    /* For rmse, use sdofminus = max(1, df_a_for_vce) when nested, df_a otherwise
       (matches ivreghdfe line 670: if HDFE.df_a=0, force absorb_ct to 1) */
    ST_int sdofminus = (has_cluster && df_a_nested > 0) ?
                       (df_a_for_vce > 0 ? df_a_for_vce : 1) : df_a;
    ST_double rmse_den = (ST_double)N_eff - K_total - n_abs_x - sdofminus
                         - sdof_extra - dofminus_vce;
    if (!small) rmse_den = N_eff - dofminus;
    ST_double rmse = sqrt(rss / (rmse_den > 0 ? rmse_den : 1));

    /* End stats timing, start store timing */
    double t_stats = ctools_timer_seconds() - t_stats_start;
    double t_store_start = ctools_timer_seconds();

    const civreghdfe_result result = {
        .df_r_val = df_r_val,
        .df_a = df_a,
        .df_a_for_vce = df_a_for_vce,
        .K_total = K_total,
        .n_abs_x = n_abs_x,
        .has_cluster2 = has_cluster2,
        .num_clusters = num_clusters,
        .num_clusters2 = num_clusters2,
        .has_cluster = has_cluster,
        .est_method = est_method,
        .K_exog_orig = K_exog_orig,
        .K_endog_orig = K_endog_orig,
        .K_endog = K_endog,
        .N_eff = N_eff,
        .rss = rss,
        .tss = tss,
        .r2 = r2,
        .rmse = rmse,
        .lambda = lambda,
        .beta = beta,
        .V = V,
        .first_stage_F = first_stage_F,
        .is_collinear_x = is_collinear_x,
        .residual_rc = residual_rc,
    };
    iv_rc = civreghdfe_post_results(&result);
    if (iv_rc) goto iv_cleanup;

    /* vce(bootstrap): keep the model for the replicates (plan 9.3) */
    if (bs_retain) {
        const ST_int K_b = K_exog_orig + K_endog_orig;
        bs_b = (ST_double *)calloc((size_t)(K_b > 0 ? K_b : 1), sizeof(ST_double));
        if (!bs_b) {
            iv_rc = 920;
            goto iv_cleanup;
        }
        /* e(b)'s layout: [endogenous, exogenous], zeros where omitted */
        for (ST_int i = 0, ce = 0, cx = 0; i < K_exog_orig + K_endog_orig; i++) {
            if (i < K_exog_orig) {
                if (!is_collinear_x[i]) bs_b[K_endog_orig + i] = beta[cx++];
            } else if (!is_collinear_x[i]) {
                bs_b[i - K_exog_orig] = beta[K_exog + ce++];
            }
        }
        const civreghdfe_boot_source bs_src = {
            .N = N, .G = G,
            .K_endog = input.K_endog, .K_exog = input.K_exog, .K_excl = K_excl_in,
            .y = y_c, .X_endog = X_endog_c, .X_exog = X_exog_c,
            .Z_excl = Z_c + (size_t)N * input.K_exog,
            .tilde = bs_tilde,
            .weights = has_weights ? weights_c : NULL, .weight_type = weight_type,
            .factors = state->factors,
            .partial = bs_partial,
            .collin_x = is_collinear_x, .collin_z = bs_collin_z, .endo_in_z = bs_endo_in_z,
            .K_b = K_b, .b_hat = bs_b,
            .est_method = est_method, .kclass_user = kclass_user, .fuller_alpha = fuller_alpha,
            .center = center, .dofminus = dofminus, .maxiter = maxiter,
            .keep_singletons = keep_singletons == 1, .tolerance = tolerance,
            .vce_type = vce_type, .kernel_type = kernel_type, .kiefer = kiefer,
            .cluster_ids = (vce_type == 2 && has_cluster) ? cluster_ids_c : NULL,
            .obs_map = est_obs,
            .N_eff = N_eff, .K_iv = K_iv,
        };
        if (!bs_tilde) iv_rc = 198;
        else iv_rc = civreghdfe_boot_retain_hook(&bs_src);
        if (iv_rc) goto iv_cleanup;
    }

    /* End store timing, compute total */
    t_store = ctools_timer_seconds() - t_store_start;
    double t_total = ctools_timer_seconds() - t_total_start;

    /* Save timing scalars (checked, as every other result) */
    {
        const char *const time_names[] = {
            "_civreghdfe_time_load", "_civreghdfe_time_extract", "_civreghdfe_time_missing",
            "_civreghdfe_time_singleton", "_civreghdfe_time_remap", "_civreghdfe_time_fwl",
            "_civreghdfe_time_partial", "_civreghdfe_time_dof", "_civreghdfe_time_estimate",
            "_civreghdfe_time_stats", "_civreghdfe_time_store", "_civreghdfe_time_total"};
        const double time_values[] = {t_load, t_extract, t_missing, t_singleton,
                                      t_hdfe_setup, t_fwl, t_partial_out, t_postproc,
                                      t_estimate, t_stats, t_store, t_total};
        for (size_t t = 0; t < sizeof(time_values) / sizeof(time_values[0]); t++) {
            if ((iv_rc = ctools_scal_save(time_names[t], time_values[t])))
                goto iv_cleanup;
        }
    }
    CTOOLS_SAVE_THREAD_INFO("_civreghdfe");

iv_cleanup:
    /* Single exit: releases whatever the stage reached had allocated,
       including after a failed result publication. The memory is returned
       by a helper thread after the plugin call (ctools_release.h): at 10M
       rows that is over 1 GB, 50 ms or more on Windows. */
    /* the tasks that return the loaded variables, before what is left of them */
    if (free_pool) {
        ctools_persistent_pool_wait(free_pool);
        free_pool = NULL;
    }
    if (prepared) ctools_release_later_call(civ_release_prepared, prepared);
    ctools_release_later(part_cols);
    ctools_release_later(cluster_ids); ctools_release_later(cluster2_ids);
    ctools_release_later(valid_mask); ctools_release_later(singleton_mask);
    if (fe_levels_c) {
        for (ST_int g = 0; g < G; g++) ctools_release_later(fe_levels_c[g]);
        ctools_release_later(fe_levels_c);
    }
    for (ST_int g = 0; g < 10; g++) {
        ctools_release_later(fe_order[g].rows.first_row);
        ctools_release_later(fe_order[g].rows.single_rows);
    }
    ctools_release_later(fe_sing.keep); ctools_release_later(fe_sing.dropped);
    for (ST_int g = 0; fe_sing.affected && g < fe_sing.G; g++) {
        ctools_release_later(fe_sing.affected[g].level);
        ctools_release_later(fe_sing.affected[g].left);
    }
    ctools_release_later(fe_sing.affected);
    for (ST_int g = 0; fe_sing.rows_left && g < fe_sing.G; g++) ctools_release_later(fe_sing.rows_left[g]);
    ctools_release_later(fe_sing.rows_left);
    ctools_release_later(level_max); ctools_release_later(is_partial);
    ctools_release_later(is_collinear_x); ctools_release_later(test_idx); ctools_release_later(ts_int);
    ctools_release_later(fe_redundant); ctools_release_later(fe_nested); ctools_release_later(fe_exact);
    for (ST_int g = 0; g < 10; g++) ctools_release_later(slope_c[g]);
    ctools_release_later(bs_partial); ctools_release_later(bs_collin_z);
    ctools_release_later(bs_tilde); ctools_release_later(bs_b);
    ctools_release_later(all_data); ctools_release_later(y_c); ctools_release_later(X_endog_c);
    ctools_release_later(X_exog_c); ctools_release_later(Z_c);
    ctools_release_later(weights_c); ctools_release_later(cluster_ids_c); ctools_release_later(cluster2_ids_c);
    ctools_release_later(beta); ctools_release_later(V);
    ctools_release_later(first_stage_F); ctools_release_later(abw_fe);
    if (state) ctools_release_later_call(civ_release_state, state);
    iv_state = NULL;
    ctools_release_later_aligned(est_obs);
    ctools_filtered_data_release(&filtered);

    return iv_rc;
}

/* The full fit (civreghdfe_main, in civreghdfe_boot.c, dispatches) */
ST_retcode civreghdfe_fit(void)
{
    return do_iv_regression();
}

/*
 * Cleanup function for civreghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 */
void civreghdfe_cleanup_state(void)
{
    /* civreghdfe cleans up iv_state at the end of do_iv_regression(),
       so this is just a safety measure for interrupted execution */
    if (iv_state != NULL) {
        /* Note: Full cleanup would require knowing G and num_threads,
           which we don't have here. The main function handles full cleanup.
           This just nulls the pointer to prevent double-free issues. */
        iv_state = NULL;
    }
    /* the model retained for vce(bootstrap) */
    if (civreghdfe_boot_release_hook) civreghdfe_boot_release_hook();
}

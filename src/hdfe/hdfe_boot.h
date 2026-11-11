/*
 * hdfe_boot.h
 *
 * Pairs-bootstrap engine shared by the HDFE estimators: resampling draws,
 * replicate scheduling, accumulation, and the summary statistics that
 * Stata's bootstrap posts. Design: docs/development/BOOTSTRAP_PLAN.md
 * (sections 5 and 6; Stata's conventions in Appendix A).
 *
 * Units and strata are 0-based C indices. Nothing here calls the Stata SPI.
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_BOOT_H
#define CTOOLS_HDFE_BOOT_H

#include <stddef.h>
#include <stdint.h>

#include "stplugin.h"

/* Replicate outcomes, in the order e(boot_status) counts them. */
enum
{
    HDFE_BOOT_OK = 0,
    HDFE_BOOT_COLLINEAR = 1, /* omitted-coefficient pattern differs from the full fit */
    HDFE_BOOT_EMPTY = 2,     /* no observations left after screening */
    HDFE_BOOT_NOCONV = 3,    /* projection or IRLS did not converge */
    HDFE_BOOT_OTHER = 4,
    HDFE_BOOT_NSTATUS = 5
};

/*
 * Resampling frame of U units in S strata. The units of stratum s are
 * stratum_start[s] .. stratum_start[s + 1] - 1 (stratum_start NULL: one
 * stratum of all U units). For draws(stata) the units must be numbered in
 * bsample's order: strata sorted, and within a stratum the rows in data
 * order or the clusters in sorted order.
 *
 * unit_freq: frequency weight of each unit, for resampling the
 * observations that frequency-weighted rows stand for (draws(fast) by
 * conditional binomials, draws(stata) through hdfe_boot_map_stata's
 * expanded frame), or NULL. Arrays are owned by the caller.
 */
typedef struct
{
    ST_int U;
    ST_int S;
    const ST_int *stratum_start;
    const ST_double *unit_freq;
    uint64_t seed; /* draws(fast) master seed, taken from Stata's RNG */
} hdfe_boot_design;

/* Returns 0, 198 for an invalid frame (a unit frequency that is not a
 * nonnegative integer of at most INT_MAX / 2, or a stratum total above 2^53),
 * or 2001 when a stratum is empty. Totals above INT_MAX are valid: draws(fast)
 * draws them; draws(stata) cannot (hdfe_boot_draw_units). */
ST_retcode hdfe_boot_design_check(const hdfe_boot_design *d);

/*
 * draws(fast): multiplicities of the U units in replicate r (0-based), a
 * pure function of (seed, r): replicate r has its own generator, and strata
 * are drawn in order. Without unit_freq, each stratum of n units draws n
 * units with replacement. With unit_freq, each stratum draws
 * Multinomial(F, freq / F), F the stratum's total frequency, by
 * conditional binomials (O(n) whatever F is).
 */
void hdfe_boot_draw_fast(const hdfe_boot_design *d, ST_int r, ST_int *unit_mult);

/*
 * draws(stata): bsample's mapping from Stata's uniforms to multiplicities
 * (Appendix A of the plan). u_draw and u_tie hold the U uniforms bsample
 * takes for the draws and for the tie-breaks, in its order. cluster_mode
 * selects bsample's cluster mapping; otherwise its observation mapping.
 * Returns 0, 198 for a uniform outside [0, 1), or 920.
 */
ST_retcode hdfe_boot_map_bsample(const hdfe_boot_design *d, const ST_double *u_draw,
                                 const ST_double *u_tie, int cluster_mode, ST_int *unit_mult);

/* The units bsample draws from under draws(stata): U, or with unit_freq the
 * observations of the fweight-expanded data (the sum of the frequencies);
 * -1 when that sum exceeds INT_MAX, which draws(stata) cannot draw. */
ST_int hdfe_boot_draw_units(const hdfe_boot_design *d);

/* The same number as a double, whatever its size (posted as __<cmd>_bs_U,
 * from which the ado sizes draws(stata)'s chunks or refuses them). */
double hdfe_boot_draw_total(const hdfe_boot_design *d);

/*
 * draws(stata) for any frame. Without unit_freq, hdfe_boot_map_bsample.
 * With unit_freq (observations resampled), bsample's observation mapping
 * on the fweight-expanded data, in the order expand leaves it and bsample
 * sorts it (strata, stable): within each stratum, the first copy of each
 * unit in unit order, then each unit's further copies in unit order.
 * u_draw and u_tie hold hdfe_boot_draw_units(d) uniforms each; unit_mult
 * gets the number of a unit's observations drawn. Returns 0, 198 or 920.
 */
ST_retcode hdfe_boot_map_stata(const hdfe_boot_design *d, const ST_double *u_draw,
                               const ST_double *u_tie, int cluster_mode, ST_int *unit_mult);

/* One replicate's results, written by the estimator's fit(). */
typedef struct
{
    ST_double *b;  /* K_b coefficients (kept regressors, then _cons) */
    ST_double *fe; /* L fixed-effect levels, NaN where not estimable; or NULL */
    ST_double *d;  /* Nd summed fixed effects per row, NaN where not identified; or NULL */
    ST_int status; /* HDFE_BOOT_* */
    ST_int iterations;
    ST_int fe_split; /* nonzero when a connected set of levels split (plan 5.7) */
    ST_int threads;  /* set by the engine: threads the fit may use for its own
                        column parallelism (OpenMP regions it opens) */
} hdfe_boot_rep;

/*
 * The estimator's side. fit() runs on a worker thread: it may not call the
 * SPI, print, or write globals. It returns 0, or a Stata error code (920)
 * that aborts the whole run; a replicate that fails statistically is
 * reported through out->status instead. Its result may not depend on
 * out->threads (column-parallel work only, each column on one thread).
 *
 * thread_bytes (optional, NULL: none): the fit's extra memory for each
 * thread beyond the first, so that the engine hands spare threads to a
 * replicate only within the budget.
 */
typedef struct
{
    size_t (*workspace_bytes)(const void *model);
    void *(*workspace_new)(const void *model);
    void (*workspace_free)(void *ws);
    ST_retcode (*fit)(const void *model, void *ws, const ST_int *unit_mult, hdfe_boot_rep *out);
    size_t (*thread_bytes)(const void *model);
} hdfe_boot_model_ops;

/*
 * Accumulated replicates. b holds every replicate's coefficients (R_cap x
 * K_b, row-major), so the summary can sort them; the fixed effects and d
 * are summarized on the fly (Welford, in replicate order).
 */
typedef struct
{
    ST_int K_b, L, Nd, R_cap;
    ST_int R; /* replicates folded so far */
    ST_double *b;
    ST_int *status;
    ST_int status_count[HDFE_BOOT_NSTATUS];
    ST_double *fe_mean, *fe_m2;
    ST_int *fe_n;
    ST_double *d_mean, *d_m2;
    ST_int *d_n;
    ST_int iterations; /* summed over replicates */
    ST_int fe_split_reps; /* replicates in which a connected set split */
    float *fe_draws;      /* L x R_cap (level-major), NaN where not estimable,
                             when kept for percentile bounds; or NULL */
} hdfe_boot_acc;

/* Returns 0 or 920; free with hdfe_boot_acc_free (idempotent). */
ST_retcode hdfe_boot_acc_init(hdfe_boot_acc *acc, ST_int K_b, ST_int L, ST_int Nd,
                              ST_int R_cap);
void hdfe_boot_acc_free(hdfe_boot_acc *acc);

/* Keep every replicate's fixed effects (in float) for percentile bounds,
 * before the first replicate is folded. Returns 0, 198 (no fixed effects
 * or replicates already folded), 909 when L x R_cap floats exceed
 * max_bytes, or 920. */
ST_retcode hdfe_boot_acc_keep_fe(hdfe_boot_acc *acc, size_t max_bytes);

/*
 * Run replicates first .. first + count - 1 and fold them into acc in
 * replicate order (acc->R must equal first). supplied: count x U unit
 * multiplicities from the ado (draws(stata), _bsdraws()), or NULL for
 * draws(fast). P = hdfe_boot_parallel(...) replicates run at once (at most
 * count), each with its own workspace, taken dynamically from the chunk
 * (from windows of a multiple of P replicates when the fixed-effect
 * outputs are large); the threads left over go to each replicate as
 * rep->threads. Results do not depend on threads. Returns 0, 1 (cancelled), 198, 909
 * (one replicate does not fit the budget), or the first error fit()
 * returned.
 */
ST_retcode hdfe_boot_run(const hdfe_boot_design *d, const ST_int *supplied,
                         const hdfe_boot_model_ops *ops, const void *model, ST_int first,
                         ST_int count, ST_int threads, size_t budget, hdfe_boot_acc *acc);

/*
 * Replicates hdfe_boot_run keeps in flight for a large chunk: threads,
 * fewer when budget bytes do not hold that many workspaces (K_b, L and Nd
 * as for hdfe_boot_acc_init); 0 when not even one fits. The ado sizes its
 * chunks as multiples of it, so that no thread idles at a chunk's end.
 */
ST_int hdfe_boot_parallel(const hdfe_boot_design *d, const hdfe_boot_model_ops *ops,
                          const void *model, ST_int K_b, ST_int L, ST_int Nd, ST_int threads,
                          size_t budget);

/* Stata's bootstrap results for the coefficients (2 x K_b intervals:
 * row 0 lower, row 1 upper; NaN where Stata posts missing). */
typedef struct
{
    ST_int K_b;
    ST_int N_reps;    /* complete replicates */
    ST_int N_misreps; /* replicates with missing coefficients */
    ST_double *b_bs, *se, *bias, *z0, *V;
    ST_double *ci_normal, *ci_percentile, *ci_bc;
} hdfe_boot_summary;

/*
 * Summarize as Stata's _bs_sum does: complete replicates only (casewise);
 * V the covariance about the replicate mean, or with mse the mean of
 * (b^r - b_hat)(b^r - b_hat)'; percentiles by _pctile's definition; the
 * bias-corrected z0 from #{b^r <= b_hat}. level in percent. Returns 0, 920,
 * or 2000 with fewer than two complete replicates.
 */
ST_retcode hdfe_boot_summarize(const hdfe_boot_acc *acc, const ST_double *b_hat, int mse,
                               ST_double level, hdfe_boot_summary *out);
void hdfe_boot_summary_free(hdfe_boot_summary *s);

/* Standard deviation over estimable replicates of each fixed-effect level
 * (or each row's d): n >= 2 gives sqrt(M2 / (n - 1)), otherwise NaN. */
void hdfe_boot_fe_se(const hdfe_boot_acc *acc, ST_double *fe_se, ST_double *d_se);

/* Percentile bounds of each fixed-effect level over the replicates that
 * estimate it, as e(ci_percentile) is computed (_pctile at (100 - level) / 2
 * and (100 + level) / 2); NaN with fewer than two. Needs the kept draws
 * (hdfe_boot_acc_keep_fe). Returns 0, 198 or 920. */
ST_retcode hdfe_boot_fe_ci(const hdfe_boot_acc *acc, ST_double level, ST_double *lb,
                           ST_double *ub);

/* _pctile's default percentile of n sorted values, p in percent. */
double hdfe_boot_pctile(const double *sorted, ST_int n, double p);

/* Standard normal CDF and its inverse (double precision). */
double hdfe_boot_normal(double x);
double hdfe_boot_invnormal(double p);

#endif /* CTOOLS_HDFE_BOOT_H */

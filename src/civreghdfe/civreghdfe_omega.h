/*
    civreghdfe_omega.h
    Covariance matrix of the IV moment conditions (ivreg2's m_omega)

    One builder serves the VCE, the efficient GMM / CUE weighting matrix, the
    Hansen J and C statistics and the Kleibergen-Paap rank statistics, so all
    of them use the same heteroskedasticity, cluster and kernel structure.
*/

#ifndef CIVREGHDFE_OMEGA_H
#define CIVREGHDFE_OMEGA_H

#include "../stplugin.h"

/* Structure of the moment covariance (the blocks of ivreg2's m_omega) */
#define CIV_OMEGA_IID     0   /* homoskedastic; with a kernel: AC */
#define CIV_OMEGA_ROBUST  1   /* heteroskedastic; with a kernel: HAC */
#define CIV_OMEGA_CLUSTER 2   /* one-way cluster */
#define CIV_OMEGA_TWOWAY  3   /* two-way cluster (Cameron-Gelbach-Miller) */
#define CIV_OMEGA_DKRAAY  4   /* clusters = time periods, kernel over periods */

/*
    tsset structure of the estimation sample, needed by the kernel
    estimators: 1-based panel ids (NULL = a single time series) and the
    0-based time index (t - tmin)/tdelta of every observation.
*/
typedef struct {
    const ST_int *panel;
    ST_int num_panels;
    const ST_int *tidx;
    ST_int tmax_idx;       /* largest time index */
    ST_double tdelta;      /* tsset delta, for the spectral-window lag count */
} civreghdfe_ts;

/*
    Builder state. civ_omega_init fills it for one estimation sample; a copy
    made with civ_omega_with_kernel shares the owned arrays (only the
    original may be passed to civ_omega_free).
*/
typedef struct {
    ST_int N;
    ST_double N_eff;            /* ivreg2's N: sum of fweights, else #obs */
    const ST_double *w;         /* weights (aw/pw normalized to sum N), NULL = none */
    ST_int wtype;               /* 0 none, 1 aw, 2 fw, 3 pw */
    ST_int kind;                /* CIV_OMEGA_* */
    ST_int kernel;              /* CIVREGHDFE_KERNEL_*, 0 = none */
    ST_double bw;
    ST_int max_lag;             /* lags 1..max_lag enter a kernel sum */
    const ST_int *cl1; ST_int ncl1;   /* 1-based cluster ids */
    const ST_int *cl2; ST_int ncl2;
    ST_int *cl12; ST_int ncl12;       /* intersection clusters (owned) */
    const ST_int *tidx;         /* time index (kernels, Driscoll-Kraay periods) */
    ST_int tmax_idx;
    ST_double tdelta;
    ST_int *order;              /* observations sorted by panel, then time (owned) */
    ST_int *pstart;             /* panel boundaries in order, npanels + 1 (owned) */
    ST_int npanels;
    ST_int center;              /* ivreg2's center: centered moments (set after init) */
    ST_int psd;                 /* ivreg2's psd0 (1) / psda (2): negative eigenvalues of S
                                   set to 0 / to their absolute value (set after init) */
    /* Bootstrap replicates (set after init; NULL/0 otherwise): mult holds
       each row's draw count m, with w = m times the fit's weight. A row's
       copies are separate observations, so the heteroskedastic terms use
       w^2 / m (aw/pw; fw and none are already linear in m). With
       mult_groups the clusters are the resampled units, each copy its own
       cluster: a cluster's outer product is divided by its count. */
    const ST_double *mult;
    ST_int mult_groups;
    ST_int owns;                /* 1 if the owned arrays belong to this copy */
} civ_omega;

/*
    Set up the builder. kind/kernel/bw select the structure; cl1/cl2 are the
    cluster ids (one-way: cl1; two-way: both; Driscoll-Kraay: none, periods
    come from ts); ts is required for kernels and Driscoll-Kraay.
    Returns 0, 920 (memory) or 198 (missing time structure).
*/
ST_retcode civ_omega_init(civ_omega *om, ST_int kind, ST_int kernel, ST_double bw,
                          ST_int N, ST_double N_eff, const ST_double *w, ST_int wtype,
                          const ST_int *cl1, ST_int ncl1,
                          const ST_int *cl2, ST_int ncl2,
                          const civreghdfe_ts *ts);

/*
    Copy of om with another kernel (0 = none), e.g. ranktest's Bartlett.
    The copy shares om's arrays and must not be freed.
*/
civ_omega civ_omega_with_kernel(const civ_omega *om, ST_int kernel);

void civ_omega_free(civ_omega *om);

/*
    S = sum of the moment contributions u_i = e_i (x) z_i with the selected
    structure; E is N x K (residual columns), Z is N x L, both column-major.
    S is (K*L) x (K*L), index k*L + l, and "raw": N times ivreg2's m_omega
    (dofminus = 0). With om->psd, S is made positive semidefinite as
    ivreg2's psd0/psda make it. Returns 0 or 920.
*/
ST_retcode civ_omega_build(const civ_omega *om, const ST_double *E, ST_int K,
                           const ST_double *Z, ST_int L, ST_double *S);

/*
    Solve the symmetric system A x = b (A n x n) with Cholesky, falling back
    to Stata's invsym generalized inverse when A is singular or a column is
    nearly dependent on the others (invsym's 1e-9 tolerance); x may alias b.
    Returns 0, or -1 when memory runs out (any symmetric A is usable).
*/
ST_int civ_sym_solve(const ST_double *A, ST_int n, const ST_double *b, ST_int nrhs, ST_double *x);

/* Eigen-decomposition of the symmetric n x n A (column-major) by cyclic
   Jacobi sweeps: eigenvalues in ev (unsorted), eigenvectors in the columns
   of Q (n x n). Returns 0, or -1 on allocation failure. */
ST_int civ_sym_eigvec(const ST_double *A, ST_int n, ST_double *ev, ST_double *Q);

/* (A)^-1 of a symmetric matrix with the same fallback. Returns 0 or -1. */
ST_int civ_sym_inverse(const ST_double *A, ST_int n, ST_double *inv);

/* Rank of scale * A (A symmetric) as Stata counts it: n minus the columns
   that invsym() zeroes (ivreg2's rankS). Returns -1 on allocation failure. */
ST_int civ_sym_rank(const ST_double *A, ST_int n, ST_double scale);

/*
    Newey-West (1994) automatic bandwidth (ivreg2's bw(auto)) for the
    Bartlett, Parzen and Quadratic Spectral kernels: e are the first-step
    residuals, zsum the row sums of the instruments ivreg2 uses, tobs the
    time span T. Returns the bandwidth (lag + 1), or -1 on failure.
*/
ST_double civ_auto_bandwidth(const ST_double *e, const ST_double *zsum, ST_int N,
                             const civreghdfe_ts *ts, ST_double tobs, ST_int kernel);

#endif /* CIVREGHDFE_OMEGA_H */

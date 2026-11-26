#ifndef CPPMLHDFE_INTERNAL_H
#define CPPMLHDFE_INTERNAL_H
#include "cppmlhdfe_options.h"
#include "cppmlhdfe_irls.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_coef.h"
#include "ctools_types.h"
#include "ctools_ols.h"
#include "ctools_config.h"
#include <float.h>
#include <limits.h>
/* Private numerical phases. Column-major sample arrays and HDFE workspaces
 * belong to the caller. Options are immutable and never read from Stata here. */
double ppml_quad_dot(const double *x, const double *y, const double *w, int n);
void ppml_uncenter_vce(double *V, const double *means, int K);
int ppml_qrsolve(const double *X, const double *y, const double *w, int N, int K, double *b,
                 double *resid);
int ppml_joint_slope(const HDFE_State *S, int g);
void ppml_project(HDFE_State *S, double *y, int g);
void ppml_transform(HDFE_State *S, const double *y, double *ans);
HDFE_SolveResult ppml_reference_partial(HDFE_State *S, double *data, int N, int K, int threads);
double ppml_norm(const double *x, int n);
void ppml_fe_multiply(HDFE_State *S, const int *offset, const double *x, double *out,
                      int transpose);
HDFE_SolveResult ppml_lsmr(const PPML_Options *options, HDFE_State *S, double *data, int N, int K);
void ppml_transform_option(const PPML_Options *options, HDFE_State *S, const double *y, double *out,
                           double *work);
HDFE_SolveResult ppml_map_options(const PPML_Options *options, HDFE_State *S, double *data, int N,
                                  int K);
HDFE_SolveResult ppml_partial(const PPML_Options *options, HDFE_State *S, double *data, int N,
                              int K, int threads);

/* What the projections of one fit keep between them while the factors stay
 * as they are and only the IRLS weights change (the pre-pass, every
 * iteration, the final projection and the refits of eta).
 *
 * plan: the coefficient solver's structure (hdfe_coef.h), built by the first
 * projection that solves for the coefficients; NULL keeps nothing.
 *
 * counts_stale: the solve for the coefficients never reads the factors'
 * weighted level counts, so the loop does not sum them after every weight
 * update: it sets this flag, and ppml_partial_kept() brings the counts up to
 * date (ppml_fe_counts_sync) before any projection that works on the rows.
 * Whoever sets the flag syncs before code that does not go through
 * ppml_partial_kept() reads the counts. */
typedef struct
{
    hdfe_coef_plan *plan;
    int counts_stale;
} ppml_projection;
/* ppml_partial() with kept (may be NULL: then exactly ppml_partial) */
HDFE_SolveResult ppml_partial_kept(const PPML_Options *options, HDFE_State *S,
                                   ppml_projection *kept, double *data, int N, int K,
                                   int threads);
/* The factors' level counts for S->weights, when kept says they lag */
void ppml_fe_counts_sync(HDFE_State *S, ppml_projection *kept);
int ppml_simplex_tableau(double *X, int n, int k, int maxiter, int *separated);
/* w_pos: the weights of the positive rows in the positive-sample
 * projection, or NULL for the fit's rule (w_user with fweights, else 1) */
int ppml_simplex(const PPML_Options *options, HDFE_State *S, const double *y, const double *X,
                 const double *w_user, int weight_type, const double *w_pos, int N, int K,
                 int *separated);
int ppml_relu_basis(HDFE_State *S, const double *X, int N, int K, double *q, int nthreads);
ST_int ppml_relu_ex(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                    const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                    ST_double *certificate);
/* ppml_relu_ex, also returning the number of ReLU iterations run (the one
 * that found the separated observations counted) in *iterations */
ST_int ppml_relu_ex2(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                     const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                     ST_double *certificate, ST_int *iterations);
/* ppmlhdfe's messages (Mata %g of a small positive number: d.dddde-XX) */
void ppml_format_small_g(char *buffer, size_t size, double value);
ST_int ppml_relu(const PPML_Options *options, HDFE_State *S, const double *y, const double *X,
                 int N, int K, int *mask, int threads);
void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N);
/* The part of ppml_update_fe_weights that the iteration on the rows reads:
 * every factor's weighted level counts (and slope centers) */
void ppml_update_fe_counts(HDFE_State *S, const ST_double *irls_w, ST_int N);

/* One IRLS update in three passes over the rows (cppmlhdfe_moments.c).
 *
 * ppml_irls_step_eta: eta_new = X_part * beta + (z - z_part), the step's
 * linear predictor without the offset; X_part (N x K) and z_part are the
 * partialled regressors and working variable. Sets *eta_not_finite when a
 * value is not finite and *mu_overflow when exp(eta_new + offset) is not
 * (neither flag is cleared).
 *
 * ppml_irls_step_means: eta = eta_new, mu = exp(eta + offset) floored at
 * PPML_MU_FLOOR, and in the same pass the largest change of eta over the
 * rows that count for convergence (the previous mu decides), whether a
 * positive outcome now has its mean at the floor, and the deviance
 * ppml_compute_deviance(y, eta, offset, w_user, N) would return (its sums,
 * in its order, from the means just computed). Rows whose eta_new is not
 * finite are left as they are and reported. Returns 0 or 920.
 *
 * ppml_irls_step_weights: irls_w = max(mu * w_user / user_weight_scale,
 * DBL_MIN), zero where excluded (may be NULL) is set, and *sum their sum over
 * the fixed blocks of ppml_update_fe_weights. Returns 0 or 920. */
typedef struct
{
    ST_double max_eta_change, deviance;
    int eta_not_finite, mu_overflow, floored_positive;
} ppml_step_means;
void ppml_irls_step_eta(const ST_double *X_part, ST_int N, ST_int K, const ST_double *beta,
                        const ST_double *z, const ST_double *z_part, const ST_double *offset,
                        ST_double *eta_new, int *eta_not_finite, int *mu_overflow);
ST_retcode ppml_irls_step_means(const ST_double *y, const ST_double *eta_new,
                                const ST_double *offset, const ST_double *w_user, ST_int N,
                                ST_double *eta, ST_double *mu, ppml_step_means *out);
ST_retcode ppml_irls_step_weights(const ST_double *mu, const ST_double *w_user,
                                  ST_double user_weight_scale, const ST_int *excluded, ST_int N,
                                  ST_double *irls_w, ST_double *sum);
/* ppmlhdfe floors the IRLS mean at Mata's epsilon(100), the spacing of
   doubles at 100 (2^-46), and has no upper bound. */
#define PPML_MU_FLOOR 1.4210854715202004e-14
/* Deviance and log likelihood from the linear predictor eta (+ offset)
   itself, mu = exp(eta + offset) without the IRLS floor; offset may be
   NULL and log_scale undoes a response scaling. */
ST_double ppml_compute_deviance(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                const ST_double *w_user, ST_int N);
ST_double ppml_compute_loglik(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              ST_double log_scale, const ST_double *w_user, ST_int N);
/* ppml_compute_loglik for eta and for the constant eta0 (no offset or
 * scale) in one pass: the same two results as two calls */
void ppml_compute_loglik_pair(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              ST_double log_scale, const ST_double *w_user, ST_double eta0, ST_int N,
                              ST_double *ll, ST_double *ll0);
/* ppml_compute_loglik_pair for the outcomes y[i] * y_scale, formed row by
 * row: the two numbers that call returns once y has been multiplied by
 * y_scale in place. Every row adds the terms it adds there; what depends on
 * the outcome alone is evaluated once per distinct value in a block of rows
 * instead of once per row. */
void ppml_compute_loglik_pair_scaled(const ST_double *y, ST_double y_scale, const ST_double *eta,
                                     const ST_double *offset, ST_double log_scale,
                                     const ST_double *w_user, ST_double eta0, ST_int N,
                                     ST_double *ll, ST_double *ll0);
/* ppml_compute_deviance(y, eta, offset, w_user, N) that also stores the
 * means, mu[i] = exp(eta[i] + offset[i]) floored at PPML_MU_FLOOR: the
 * deviance takes each exp() as evaluated for mu (before the floor), so both
 * are the numbers of the two separate passes. */
ST_double ppml_means_deviance(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              const ST_double *w_user, ST_int N, ST_double *mu);
/* The double-double sum of x over fixed blocks of rows, each in row order,
 * the blocks then added in order: the same at every thread count */
dd_real ppml_block_sum_dd(const ST_double *x, ST_int N);
/* fe = eta - X b for the N x K column-major X. X b is left in xb (N): every
 * row accumulates its products over the regressors in order from zero, as
 * ctools_matvec_add() into a zeroed vector does. */
void ppml_fixed_effect_part(const ST_double *X, ST_int N, ST_int K, const ST_double *beta,
                            const ST_double *eta, ST_double *xb, ST_double *fe);
#endif

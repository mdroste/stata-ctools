/*
 * cqreg_sparsity.h
 *
 * Sparsity estimation for quantile regression variance computation.
 * Implements bandwidth selection and kernel density estimation.
 * Part of the ctools suite.
 */

#ifndef CQREG_SPARSITY_H
#define CQREG_SPARSITY_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Interface
 * ============================================================================ */

/*
 * Estimate sparsity (1/f(F^{-1}(q))) from regression residuals.
 *
 * Parameters:
 *   sp        - Pre-allocated sparsity state
 *   residuals - Regression residuals (N)
 *
 * Returns:
 *   Estimated sparsity value (stored in sp->sparsity)
 */
ST_double cqreg_estimate_sparsity(cqreg_sparsity_state *sp,
                                  const ST_double *residuals);

/* ============================================================================
 * Bandwidth Selection Methods
 * ============================================================================ */

/*
 * Hall-Sheather (1988) bandwidth.
 * Optimal for kernel density estimation at a quantile.
 *
 * h = z_{1-alpha/2}^{2/3} * (1.5 * phi(z_q)^2 / (2*z_{1-alpha/2}^2 + 1))^{1/3} * n^{-1/3}
 *
 * where z_q = Phi^{-1}(q), phi = standard normal PDF
 */
ST_double cqreg_bandwidth_hsheather(ST_int N, ST_double q, ST_double alpha);

/*
 * Bofinger (1975) bandwidth.
 *
 * h = (9/2 * phi(z_q)^4 * (2*z_q^2 + 1)^2 / n)^{1/5}
 */
ST_double cqreg_bandwidth_bofinger(ST_int N, ST_double q);

/*
 * Chamberlain (1994) bandwidth.
 *
 * h = z_{1-alpha/2} * sqrt(q*(1-q) / n) / phi(Phi^{-1}(q))
 */
ST_double cqreg_bandwidth_chamberlain(ST_int N, ST_double q, ST_double alpha);

/*
 * Compute bandwidth using the specified method.
 */
ST_double cqreg_compute_bandwidth(ST_int N, ST_double q, cqreg_bw_method method);

/* ============================================================================
 * Kernel Functions
 * ============================================================================ */

/*
 * Epanechnikov kernel: K(u) = 0.75 * (1 - u^2) for |u| < 1, 0 otherwise.
 * Optimal in MISE sense.
 */
ST_double cqreg_kernel_epanechnikov(ST_double u);

/*
 * Gaussian kernel: K(u) = (1/sqrt(2*pi)) * exp(-u^2/2)
 */
ST_double cqreg_kernel_gaussian(ST_double u);

/*
 * Triangular kernel: K(u) = (1 - |u|) for |u| < 1, 0 otherwise.
 */
ST_double cqreg_kernel_triangular(ST_double u);

/* ============================================================================
 * Density Estimation
 * ============================================================================ */

/*
 * Estimate density at a point using kernel density estimation.
 *
 * f_hat(x) = (1/(n*h)) * sum_i K((x - X_i) / h)
 *
 * Parameters:
 *   x         - Point at which to estimate density
 *   data      - Data points (N)
 *   N         - Number of data points
 *   bandwidth - Kernel bandwidth
 *
 * Returns:
 *   Estimated density at x
 */
ST_double cqreg_kernel_density(ST_double x,
                               const ST_double *data,
                               ST_int N,
                               ST_double bandwidth);

/*
 * Estimate density at the q-th quantile of the data.
 * This is what's needed for sparsity estimation.
 *
 * Parameters:
 *   data      - Data points (N)
 *   N         - Number of data points
 *   q         - Quantile (0 < q < 1)
 *   bandwidth - Kernel bandwidth
 *
 * Returns:
 *   Estimated density at the q-th quantile
 */
ST_double cqreg_density_at_quantile(const ST_double *data,
                                    ST_int N,
                                    ST_double q,
                                    ST_double bandwidth);

/* ============================================================================
 * Statistical Functions
 * ============================================================================ */

/*
 * Standard normal CDF: Phi(x)
 */
ST_double cqreg_pnorm(ST_double x);

/*
 * Standard normal PDF: phi(x)
 */
ST_double cqreg_dnorm(ST_double x);

/*
 * Standard normal quantile function: Phi^{-1}(p)
 */
ST_double cqreg_qnorm(ST_double p);

/*
 * Compute sample quantile.
 * Uses linear interpolation (type 7 in R terminology).
 *
 * Parameters:
 *   data   - Sorted data (N)
 *   N      - Number of observations
 *   q      - Quantile (0 < q < 1)
 *
 * Returns:
 *   Sample quantile value
 */
ST_double cqreg_sample_quantile(const ST_double *sorted_data, ST_int N, ST_double q);

/*
 * Sort array in ascending order (in-place).
 */
void cqreg_sort_ascending(ST_double *data, ST_int N);

/* ============================================================================
 * Alternative Density Estimation Methods
 * ============================================================================ */

/*
 * Estimate per-observation densities using the "fitted" method.
 * This is Stata's default density estimation for qreg.
 *
 * Parameters:
 *   obs_density   - Output: per-observation density estimates (N)
 *   residuals     - Regression residuals (N)
 *   N             - Number of observations
 *   q             - Quantile
 *   bw_method     - Bandwidth selection method
 *
 * Returns:
 *   Average sparsity (1 / mean density)
 */
ST_double cqreg_estimate_fitted_density(ST_double *obs_density,
                                        const ST_double *residuals,
                                        ST_int N,
                                        ST_double q,
                                        cqreg_bw_method bw_method);

/*
 * Estimate sparsity using kernel density estimation on residuals.
 */
ST_double cqreg_estimate_kernel_density_sparsity(cqreg_sparsity_state *sp,
                                                  const ST_double *residuals);

#endif /* CQREG_SPARSITY_H */

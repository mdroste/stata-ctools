#ifndef CPPMLHDFE_INTERNAL_H
#define CPPMLHDFE_INTERNAL_H
#include "cppmlhdfe_options.h"
#include "cppmlhdfe_irls.h"
#include "hdfe/hdfe_solver.h"
#include "ctools_types.h"
#include "ctools_ols.h"
#include "ctools_config.h"
#include <float.h>
#include <limits.h>
/* Private numerical phases. Column-major sample arrays and HDFE workspaces
 * belong to the caller. Options are immutable and never read from Stata here. */
double ppml_quad_dot(const double *x, const double *y, const double *w, int n);
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
int ppml_simplex_tableau(double *X, int n, int k, int maxiter, int *separated);
int ppml_simplex(const PPML_Options *options, HDFE_State *S, const double *y, const double *X,
                 const double *w_user, int weight_type, int N, int K, int *separated);
int ppml_relu_basis(HDFE_State *S, const double *X, int N, int K, double *q, int nthreads);
ST_int ppml_relu_ex(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                    const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                    ST_double *certificate);
ST_int ppml_relu(const PPML_Options *options, HDFE_State *S, const double *y, const double *X,
                 int N, int K, int *mask, int threads);
void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N);
ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu, const ST_double *w_user,
                                ST_int N);
ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu, const ST_double *w_user,
                              ST_int N);
#endif

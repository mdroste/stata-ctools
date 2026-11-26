/*
 * hdfe_coef.h
 *
 * Fixed-effect projection solved for the coefficients of the effects
 * Part of the ctools Stata plugin suite
 *
 * The row-space iteration of hdfe_solver.c runs conjugate gradients on
 * N-length vectors. The same least-squares problem can be solved on the
 * levels of the effects: the effect with the most levels is eliminated
 * exactly (it needs only its level means), and the coefficients of the
 * other effects solve a system whose size is their number of levels,
 *
 *     S b = Dk' W M_B x,     S = Dk' W M_B Dk,
 *
 * with Dk the dummies of the kept effects, W the weights and M_B the
 * annihilator of the eliminated effect. For two effects S is the weighted
 * Laplacian of the graph that links two levels of the kept effect whenever
 * a level of the eliminated effect meets both.
 *
 * S is solved
 *   - by a dense Cholesky factorization when the kept effects have few
 *     levels (year, industry, state): no iteration and no tolerance;
 *   - otherwise by conjugate gradients on the distinct cells of the effects,
 *     for all columns in one sweep. Preconditioned by the level counts this
 *     would be, in exact arithmetic, the row-space iteration; the diagonal
 *     of S (two effects) takes a quarter to a third fewer steps;
 *   - for two effects whose graph is poorly connected (workers who rarely
 *     change firm), by conjugate gradients preconditioned with an
 *     approximate Cholesky factor of the Laplacian (Kyng and Sachdeva 2016),
 *     which needs a few dozen iterations whatever the graph.
 *
 * Every sum is accumulated over fixed blocks combined in a fixed order, so
 * results do not depend on the thread count or the schedule. The module has
 * no static data and allocates its own scratch: separate states can be
 * projected concurrently, and parallel regions take their thread count from
 * the calling context (the command's threads, or a bootstrap replicate's
 * share).
 *
 * Not for weights that fall without bound, as the IRLS weights of a Poisson
 * fit do for zero outcomes whose mean heads to zero. An effect tied to the
 * others only through rows that weigh 1e-11 of its other rows has a pivot
 * below rounding: the factorization holds it at zero and the iteration never
 * moves it, where the row-space iteration still carries those rows along.
 * The weighted fit is the same to the tolerance, but the residuals of those
 * rows are not, and cppmlhdfe needs them (cppmlhdfe_solvers.c).
 * hdfe_coef_applies() therefore declines weights whose positive values span
 * more than HDFE_COEF_WEIGHT_SPAN; partial_out_columns() then iterates on
 * the rows. Weights of exactly zero take their rows out and are fine.
 */

#ifndef CTOOLS_HDFE_COEF_H
#define CTOOLS_HDFE_COEF_H

#include "hdfe_solver.h"

typedef enum
{
    /* Dense when the kept effects are few; otherwise iterate, and for two
     * effects switch to the Laplacian preconditioner when the plain
     * iteration has not converged after HDFE_COEF_PROBE_ITERATIONS. */
    HDFE_COEF_AUTO = 0,
    HDFE_COEF_DENSE,  /* dense Cholesky of S */
    HDFE_COEF_PCG,    /* conjugate gradients: diagonal of S (two effects) or
                         symmetric block Gauss-Seidel (three or more) */
    HDFE_COEF_PCG_LAP /* two effects: approximate-Cholesky preconditioner */
} hdfe_coef_method;

/* Largest number of effects solved here; more keep the row-space iteration */
#define HDFE_COEF_MAX_EFFECTS 8

/* partial_out_columns() takes an iteration limit below this as a limit on
 * the row-space iteration and runs that one. reghdfe's iterate() counts its
 * sweeps over the rows; a direct solve has none to count, so iterate(1)
 * could never stop it, and the result after iterate(2) abort(0) is defined
 * by those sweeps alone. Limits users raise for hard problems (tens of
 * thousands) stay above it. */
#define HDFE_COEF_MIN_MAXITER 1000

/* Plain iterations tried before the Laplacian preconditioner is built. Well
 * connected data converge in 3 to 50 (validation/performance/
 * benchmark_hdfe_coef.c: random, panel and regional worker-firm graphs at
 * 10M rows need 3 to 41); data that need more gain a factor of 10 or more
 * from the preconditioner, which costs about as much as these iterations. */
#define HDFE_COEF_PROBE_ITERATIONS 60

/* 1 when S can be projected here: two to HDFE_COEF_MAX_EFFECTS factors, all
 * of them intercepts (no slope). */
/* Largest ratio of the heaviest to the lightest positive weight the module
 * takes; beyond it hdfe_coef_applies() is 0 and the rows are iterated */
#define HDFE_COEF_WEIGHT_SPAN 1e8
int hdfe_coef_applies(const HDFE_State *S);

/*
 * Residuals of the (weighted) projection of K columns on S's factors.
 *
 *   data: column-major N x K, overwritten with the residuals
 *   method: HDFE_COEF_AUTO in production; the others select one path
 *
 * Uses S->factors[g].levels and num_levels, S->weights when S->has_weights,
 * S->tolerance and S->maxiter. Returns the status of partial_out_columns():
 * 0, 1 (Break), 198 (invalid state), 920 (memory) or 430 (no convergence
 * in S->maxiter iterations; data then holds the last iterate, as the
 * row-space solver leaves it). iterations is 1 for a direct solve.
 *
 * The stopping rule is that of the row-space iteration stated on the exact
 * error: a column stops when its last step reduced the squared (weighted)
 * distance to the projection by at most tolerance^2 times the residual sum
 * of squares. A column whose solve reaches the rounding of its own size
 * first (level means far above the within variation, or a tolerance near
 * 1e-13) is solved once more from its residuals, so the tolerance is met
 * relative to what is left and not to what was removed.
 */
HDFE_SolveResult hdfe_coef_partial_out(const HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                       hdfe_coef_method method);

#endif /* CTOOLS_HDFE_COEF_H */

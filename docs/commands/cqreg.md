# cqreg

Quantile regression.

## Overview

`cqreg` estimates quantile regressions using an interior-point solver written in C. See the [published benchmarks](../benchmarks/readme_profile/results.md) for timing comparisons with `qreg`. `absorb()` is currently unsupported; use explicit factor-variable indicators to estimate fixed effects jointly with the slopes.

`fweight`s, `iweight`s and `pweight`s are accepted as in `qreg` (see [Weights](#weights)). Solver and VCE differences are described in the [compatibility table](../COMPATIBILITY.md).

## Syntax

```stata
cqreg depvar [indepvars] [if] [in] [weight] [, options]
cqreg [, level(#) display_options]
```

Without `indepvars` the constant-only model is fitted. Typed without a variable list, or through
`estimates replay`, `cqreg` redisplays the last estimates.

## Options

### Model Options
| Option | Description |
|--------|-------------|
| `quantile(#)` | Quantile to estimate (0-1, default: 0.5 for median) |
| `absorb(varlist)` | Currently unsupported; use explicit indicators |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce([vcetype] [, vceopts])` | Variance estimation method: `iid` (default; `robust` with `pweight`s), `robust` (or `r`), or `cluster clustvar` (or `qreg`'s `excluster clustvar`). `vceopts` are `qreg`'s density and bandwidth suboptions, e.g. `vce(robust, kernel(gaussian) bofinger)`, `vce(iid, res)`, `vce(, residual)` |
| `denmethod(method)` | Density estimation: `fitted` (default), `residual` or `residuals` (with `vce(iid)` only), or `kernel[(kernel)]`, with `qreg`'s abbreviations (`fit`, `res`, `ker()`); the same as the `vce()` suboptions, with which it must agree |
| `bwmethod(method)` | Bandwidth selection: `hsheather` (default), `bofinger`, or `chamberlain`, or `qreg`'s `hs`, `bo`, `ch` |

### Optimization Options
| Option | Description |
|--------|-------------|
| `tolerance(#)` | Convergence tolerance (default: 1e-12) |
| `maxiter(#)` | Maximum IPM iterations (default: 200) |
| `nopreprocess(#)` | Accepted for backward compatibility; has no effect (`nopreprocess(-1)` returns error 198) |
| `wlsiter(#)`, `iterate(#)`, `trace` | `qreg`'s options for its simplex solver (WLS starting iterations, its iteration limit, a coefficient trace), checked as `qreg` checks them (error 198 unless # is a positive integer). The interior-point solver uses neither starting values nor `iterate()`; `trace` shows the log line and the solution's coefficients |

Under `version 12` or earlier, `qreg` reports its old standard errors; `cqreg` reports those of
`qreg` 13 and later and says so in a note.

If the solver reaches its iteration limit before converging, encounters a numerical failure, or produces nonfinite results, the command returns error 430 and does not store new estimates.

### Reporting Options
| Option | Description |
|--------|-------------|
| `level(#)` | Confidence level (default 95); as in `qreg`, alpha = 1 - #/100 also enters the Hall-Sheather and Chamberlain bandwidths |
| display options | `qreg`'s estimation display options, such as `cformat()` |
| `nolog` | Hides the iteration log line, `Iteration N:  Sum of abs. weighted deviations = ...` (the interior-point iteration count and the minimized sum, as `qreg`'s last log line); `log` shows it, and `set iterlog off` hides it as in `qreg` |
| `verbose` | Display progress information |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |

## Examples

```stata
* Basic median regression (default quantile = 0.5)
sysuse auto, clear
cqreg price mpg weight

* 25th percentile regression
cqreg price mpg weight, quantile(0.25)

* 75th percentile with robust standard errors
cqreg price mpg weight, quantile(0.75) vce(robust)

* 90th percentile with clustered standard errors
cqreg price mpg weight, quantile(0.90) vce(cluster foreign)

* Quantile regression with fixed effects
webuse nlswork, clear
cqreg ln_wage age ttl_exp tenure i.idcode

* Two-way fixed effects with clustering
cqreg ln_wage age ttl_exp tenure i.idcode i.year, vce(cluster idcode)

* Using Bofinger bandwidth
sysuse auto, clear
cqreg price mpg weight, quantile(0.5) bwmethod(bofinger)

* With verbose timing
cqreg price mpg weight, verbose
```

## Performance

`cqreg` achieves its speed through:

### Implementation details

**IPM Solver:**
- **Mehrotra predictor-corrector**: Each Newton step computes a predictor and a corrector direction. Runtime depends on sample size, the number of regressors, conditioning, and the iterations needed to converge.

**Linear Algebra:**
- **BLAS/LAPACK dispatch**: Automatically uses Apple Accelerate (macOS) or OpenBLAS when available for matrix-matrix and matrix-vector operations, and uses C routines when those libraries are unavailable
- **K-way unrolled dot product**: `ctools_dot_unrolled` (8/16-way) used in the fallback path when BLAS is unavailable
- **4-way unrolled scaling**: Custom `dscal`-equivalent loop processes 4 elements per iteration when BLAS is not linked
- **OpenMP parallelization**: Parallel matrix operations in the IPM core

## Algorithm

Quantile regression minimizes:
```
minimize Σᵢ [ q * uᵢ + (1-q) * vᵢ ]
subject to: yᵢ - xᵢ'β = uᵢ - vᵢ, uᵢ ≥ 0, vᵢ ≥ 0
```

where:
- `q` is the quantile (0 < q < 1)
- `uᵢ` represents positive deviations
- `vᵢ` represents negative deviations

The IPM solver converts this to a barrier problem and uses Newton steps with predictor-corrector improvements.

## Standard Errors

### IID (Default)
Assumes independent and identically distributed errors:
```
V = sparsity² * q*(1-q) * (X'X)⁻¹
```

### Robust
Uses the Powell sandwich estimator for heteroskedasticity-robust inference. As in `qreg`, the
Hessian is inverted with `invsym()`, so directions without density (a level or the constant whose
observations all have zero fitted density) get zero rows in `e(V)`, and `e(df_r)` is N - `e(rank)`.

### Density Methods
The sparsity 1/f(0) (or, for robust standard errors, the per-observation densities) is estimated as in `qreg`:

| Method | Description |
|--------|-------------|
| `fitted` | Default. Quantile regressions at q-h and q+h, evaluated at the mean of X (iid) or at each observation (robust) |
| `residual` | Difference quotient of residual percentiles at q±h after dropping the K observations of the linear programming basis; `vce(iid)` only (error 184 otherwise) |
| `kernel[(kernel)]` | Kernel estimate from the residuals with bandwidth min(sd, IQR/1.34)·(Φ⁻¹(q+h) − Φ⁻¹(q−h)); kernels `epanechnikov` (default), `epan2`, `biweight`, `cosine`, `gaussian`, `parzen`, `rectangle`, `triangle` |

### Weights
As in `qreg`, `fweight`s, `iweight`s and `pweight`s are allowed, and the coefficients minimize
the weighted check loss Σᵢ wᵢ ρ_q(rᵢ). Because ρ_q(w·r) = w·ρ_q(r) for w > 0, the solver fits
the rows of y and X scaled by the weights. `e(sum_adev)`, `e(sum_rdev)` and `e(q_v)` are
weighted (`e(q_v)` is the first value in sort order whose cumulative weight reaches q times the
total), `e(sum_w)` is the sum of the weights, and `e(wtype)` and `e(wexp)` are stored.

- `fweight`: `e(N)`, `e(df_r)` and the bandwidth use the sum of the weights; X'X is
  Σ wᵢ xᵢxᵢ'. The residual method removes one copy of each basis observation.
- `iweight`: `e(N)` and the bandwidth count observations; X'X is Σ wᵢ² xᵢxᵢ'. The residual
  method drops the basis observations.
- `pweight`: fitted as `iweight`s; imply `vce(robust)`. `vce(iid)` and `denmethod(residual)`
  return error 404, as `qreg`'s `vce(iid)` and `vce(, residual)`.

The robust Hessian weights each fitted density or kernel value by wᵢ, the fitted iid sparsity
uses the weighted mean of X, the kernel method uses the weighted mean, standard deviation and
interquartile range of the residuals, and the cluster scores are weighted by wᵢ. Observations
with zero weight are excluded. `qreg` also accepts negative `iweight`s, for which the weighted
check loss is not convex and has no well-defined minimum; `cqreg` rejects them with error 402.

### Clustered
Cluster-robust standard errors clustered on the specified variable. With a single cluster the
covariance is undefined; as `qreg`, `cqreg` then posts a zero `e(V)` with `e(rank)` 0, so
`e(df_r)` is N and no standard errors are displayed.

### Zero sparsity
When the dependent variable has a mass point at the fitted quantile, both auxiliary fits of the
fitted method give the same quantile at the mean of X. `qreg`'s exact vertices then give a zero
sparsity and a zero covariance. `cqreg`'s interior-point auxiliary fits leave a rounding-level
difference, so a difference below sqrt(epsilon) in the scaled response units is treated as zero:
the cutoff `qreg` applies to the per-observation robust densities, which `cqreg` also applies in
those units.

### Non-unique optima
`cqreg` checks whether coefficients can change without increasing the
quantile-regression objective. If it finds such a change, it reports a nonunique
solution and sets `e(nonunique)` to 1. The check does not detect every nonunique
solution, so a value of 0 does not establish uniqueness.

### Bandwidth Methods
| Method | Description |
|--------|-------------|
| `hsheather` | Hall–Sheather (1988); default |
| `bofinger` | Bofinger (1975) |
| `chamberlain` | Chamberlain (1994) |

## Stored Results

`cqreg` stores the following in `e()`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(q)` | Quantile estimated |
| `e(sum_adev)` | Sum of absolute deviations |
| `e(sparsity)` | Estimated sparsity (1/f(0)); missing with a robust or cluster VCE, as in `qreg` |
| `e(bwidth)` | Bandwidth used for sparsity estimation |
| `e(kbwidth)` | Kernel bandwidth (kernel method) |
| `e(iterations)` | Number of IPM iterations |
| `e(convcode)` | 0 for a successfully converged fit |
| `e(df_r)` | Residual degrees of freedom, N - `e(rank)` |
| `e(rank)` | Rank of `e(V)`, assessed in the scaled response units |
| `e(df_m)` | Model degrees of freedom |
| `e(N_clust)` | Number of clusters (if clustered) |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `cqreg` |
| `e(depvar)` | Dependent variable name |
| `e(vce)` | Variance estimation method |
| `e(vcetype)` | `Robust` with a robust or cluster VCE |
| `e(denmethod)` | Density estimation method |
| `e(kernel)` | Kernel function (kernel method) |
| `e(bwmethod)` | Bandwidth selection method |
| `e(clustvar)` | Cluster variable (if clustered) |
| `e(predict)` | Program used for `predict` |
| `e(basis)` | Hidden, as in `qreg`: observation numbers of the linear programming basis (empty if the solver identified none) |

`e(tolerance)` (hidden) is the solver's `tolerance()`; `qreg`'s simplex posts its own, 1e-10. As
`qreg`, `cqreg` sets the legacy globals `$S_E_mdf`, `$S_E_tdf`, `$S_E_rd`, `$S_E_nobs`, `$S_E_q`,
`$S_E_rq`, `$S_E_rsd`, `$S_E_msd`, `$S_E_vl`, `$S_E_frc` and `$S_E_cmd` (`cqreg`); `$S_E_vl` lists
the coefficients' names, where `qreg` uses `_rmcoll`'s compact form for factor variables.

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix |

## Technical Notes

### Fixed Effects
`absorb()` returns an error because least-squares partialling does not preserve
quantile loss. Use explicit indicators when the number of groups permits it.

### Scale of the regressors
The collinearity screen runs on the correlation matrix of the regressors, as `_rmcoll`, and the
diagonal ridge of the solver's normal equations is at most 1e-14 times a column's sum of squares,
so a regressor measured on a very small scale (sum of squares near 1e-13) is neither omitted nor
shrunk. `e(rank)` is computed with `invsym()` on the covariance in the scaled response units
(spread of y between 1 and 2), not on the covariance divided by its largest entry, so a constant
beside such a regressor keeps its rank. At scales below about 1e-10 `qreg` itself omits the
regressor or stops with r(2001); `cqreg` gives the rescaled fit.

### Multiple optima
The quantile regression objective can have a set of minimizers, which is common with indicator
regressors: when a group's size times q is an integer, the group's coefficient can move within an
interval at the same objective. `qreg` and `cqreg` then can report different, equally optimal
coefficients with the same `e(sum_adev)`. `validation/stata/cqreg.do` checks such models with
`benchmark_qreg_altopt`, which verifies that cqreg's coefficients attain qreg's minimum.

### Solver
`cqreg` uses the Frisch–Newton interior-point solver. `nopreprocess(0)` and `nopreprocess(1)` are accepted for compatibility and do not change the solver.

## Comparison with qreg

| Feature | Stata `qreg` | `cqreg` |
|---------|--------------|---------|
| Algorithm | Simplex | Interior Point Method |
| HDFE Support | No | No; explicit indicators required |
| Parallelization | No | Yes (OpenMP) |

## References

- Koenker, R. 2005. *Quantile Regression*. Cambridge University Press.
- Chernozhukov, V., I. Fernández-Val, and B. Melly. 2020. Fast algorithms for the quantile regression process. *Empirical Economics* 62: 7-33.
- Portnoy, S., and R. Koenker. 1997. The Gaussian hare and the Laplacian tortoise. *Statistical Science* 12: 279-300.
- Hall, P., and S. J. Sheather. 1988. On the distribution of a studentized quantile. *JRSS-B* 50: 381-391.

## See Also

- [ctools Overview](../../README.md)
- [creghdfe](creghdfe.md) - Linear regression with HDFE
- Stata's [qreg](https://www.stata.com/manuals/rqreg.pdf) - Native quantile regression

# cqreg

High-performance C-accelerated quantile regression.

## Overview

`cqreg` is a drop-in replacement for Stata's `qreg` command that uses an Interior Point Method (IPM) solver implemented in C. It provides significant performance improvements, especially for large datasets (N > 100,000). Unlike native `qreg`, `cqreg` supports high-dimensional fixed effects via the `absorb()` option.

## Syntax

```stata
cqreg depvar indepvars [if] [in] [, options]
```

## Options

### Model Options
| Option | Description |
|--------|-------------|
| `quantile(#)` | Quantile to estimate (0-1, default: 0.5 for median) |
| `absorb(varlist)` | Categorical variables to absorb as fixed effects (HDFE) |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce(vcetype)` | Variance estimation method: `iid` (default), `robust`, or `cluster clustvar` |
| `bwmethod(method)` | Bandwidth selection: `hsheather` (default), `bofinger`, or `chamberlain` |

### Optimization Options
| Option | Description |
|--------|-------------|
| `tolerance(#)` | Convergence tolerance (default: 1e-8) |
| `maxiter(#)` | Maximum IPM iterations (default: 50) |
| `nopreprocess(#)` | Preprocessing control: 0 = disabled (default), -1 = enabled (experimental) |

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information |
| `timeit` | Display execution time |

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
cqreg ln_wage age ttl_exp tenure, absorb(idcode)

* Two-way fixed effects with clustering
cqreg ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)

* Using Bofinger bandwidth
cqreg price mpg weight, quantile(0.5) bwmethod(bofinger)

* With verbose timing
cqreg price mpg weight, verbose timeit
```

## Performance

`cqreg` achieves its speed through:

### Speedup Tricks

**IPM Solver:**
- **Mehrotra predictor-corrector**: Each Newton step computes both a predictor (affine) and corrector direction, achieving superlinear convergence and typically finishing in 15-30 iterations regardless of N
- **O(sqrt(N)) convergence**: IPM complexity is polynomial, compared to the simplex method's exponential worst case—dominant advantage on large datasets

**Linear Algebra:**
- **BLAS/LAPACK dispatch**: Automatically uses Apple Accelerate (macOS) or OpenBLAS when available for matrix-matrix and matrix-vector operations, falling back to hand-tuned C kernels
- **K-way unrolled dot product**: `ctools_dot_unrolled` (8/16-way) used in the fallback path when BLAS is unavailable
- **4-way unrolled scaling**: Custom `dscal`-equivalent loop processes 4 elements per iteration when BLAS is not linked
- **OpenMP parallelization**: Parallel matrix operations in both the IPM core and the HDFE pre-processing step

**HDFE Support:**
- **CG solver with symmetric Kaczmarz**: Same algorithm as `creghdfe`—fixed effects are partialled out before the IPM solve, so the IPM operates on a smaller effective problem
- **Pre-computed inverse group counts**: Avoids division in the inner demeaning loop
- **Thread-local accumulators**: Eliminates contention during parallel FE projection

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
Assumes independent and identically distributed errors. Uses kernel density estimation for the sparsity function f(0):
```
V = (1/n) * sparsity² * q*(1-q) * (X'X)⁻¹
```

### Robust
Uses the Powell sandwich estimator for heteroskedasticity-robust inference.

### Clustered
Cluster-robust standard errors clustered on the specified variable.

### Bandwidth Methods
| Method | Description |
|--------|-------------|
| `hsheather` | Hall-Sheather (1988) - default, generally recommended |
| `bofinger` | Bofinger (1975) - more conservative |
| `chamberlain` | Chamberlain (1994) - less smoothing |

## Stored Results

`cqreg` stores the following in `e()`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(q)` | Quantile estimated |
| `e(sum_adev)` | Sum of absolute deviations |
| `e(sparsity)` | Estimated sparsity (1/f(0)) |
| `e(bwidth)` | Bandwidth used for sparsity estimation |
| `e(iterations)` | Number of IPM iterations |
| `e(converged)` | 1 if converged, 0 otherwise |
| `e(df_r)` | Residual degrees of freedom |
| `e(df_m)` | Model degrees of freedom |
| `e(df_a)` | Degrees of freedom absorbed by FEs |
| `e(N_clust)` | Number of clusters (if clustered) |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `cqreg` |
| `e(depvar)` | Dependent variable name |
| `e(vce)` | Variance estimation method |
| `e(bwmethod)` | Bandwidth selection method |
| `e(absorb)` | Absorbed FE variables |
| `e(clustvar)` | Cluster variable (if clustered) |
| `e(predict)` | Program used for `predict` |

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix |

## Technical Notes

### Fixed Effects (HDFE)
When `absorb()` is specified:
- Fixed effects are partialled out using CG solver before IPM
- Absorbed factors consume degrees of freedom
- More efficient than creating dummy variables

### Preprocessing (Experimental)
The `nopreprocess(-1)` option enables an experimental preprocessing algorithm based on Chernozhukov, Fernández-Val, and Melly (2020). This attempts to speed up estimation by:
1. Solving on a subsample first
2. Exploiting that most residual signs are predictable

**Note**: The preprocessing implementation is under development and disabled by default. The direct Frisch-Newton solver is recommended.

## Comparison with qreg

| Feature | Stata `qreg` | `cqreg` |
|---------|--------------|---------|
| Algorithm | Simplex | Interior Point Method |
| Complexity | O(N²) worst case | O(sqrt(N)) |
| HDFE Support | No | Yes |
| Parallelization | No | Yes (OpenMP) |
| Typical Speedup | 1x (baseline) | 2-10x |

## References

- Koenker, R. 2005. *Quantile Regression*. Cambridge University Press.
- Chernozhukov, V., I. Fernández-Val, and B. Melly. 2020. Fast algorithms for the quantile regression process. *Empirical Economics* 62: 7-33.
- Portnoy, S., and R. Koenker. 1997. The Gaussian hare and the Laplacian tortoise. *Statistical Science* 12: 279-300.
- Hall, P., and S. J. Sheather. 1988. On the distribution of a studentized quantile. *JRSS-B* 50: 381-391.

## See Also

- [ctools Overview](../README.md)
- [creghdfe](README_creghdfe.md) - Linear regression with HDFE
- Stata's [qreg](https://www.stata.com/manuals/rqreg.pdf) - Native quantile regression

# civreghdfe

High-performance C-accelerated instrumental variables regression with high-dimensional fixed effects.

## Overview

`civreghdfe` is a high-performance replacement for `ivreghdfe` that uses a C plugin for all computations. It performs instrumental variables (IV) regression with two-stage least squares (2SLS) while absorbing high-dimensional fixed effects. The command combines the speed of the `creghdfe` HDFE absorption algorithm with 2SLS estimation.

## Syntax

```stata
civreghdfe depvar (endogvars = instruments) [exogvars] [if] [in] [weight], absorb(varlist) [options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `absorb(varlist)` | **Required.** Categorical variables to absorb as fixed effects |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce(vcetype)` | Variance estimation: `unadjusted`, `robust`, or `cluster clustvar` |

### Estimation Options
| Option | Description |
|--------|-------------|
| `tolerance(#)` | CG solver convergence tolerance (default: 1e-8) |
| `maxiter(#)` | Maximum CG iterations (default: 500) |

### Reporting Options
| Option | Description |
|--------|-------------|
| `first` | Report first-stage regression statistics |
| `verbose` | Display progress information |
| `timeit` | Display timing breakdown |
| `small` | Use small-sample adjustments |

### Weights
Supports `aweight`, `fweight`, and `pweight`.

## Examples

```stata
* Basic IV regression with one endogenous variable
webuse nlswork, clear
civreghdfe ln_wage (tenure = union) age ttl_exp, absorb(idcode)

* Multiple endogenous variables
civreghdfe ln_wage (tenure hours = union south) age, absorb(idcode year)

* With robust standard errors
civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(robust)

* With clustered standard errors
civreghdfe ln_wage (tenure = union) age, absorb(idcode year) vce(cluster idcode)

* Display first-stage statistics
civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) first

* With weights
civreghdfe ln_wage (tenure = union) age [aw=hours], absorb(idcode)

* Multiple fixed effects
civreghdfe y (x1 = z1 z2) x2 x3, absorb(firm year industry)

* Verbose timing output
civreghdfe ln_wage (tenure = union) age, absorb(idcode) verbose timeit
```

## Syntax Details

The endogenous variables and their instruments are specified in parentheses:
```
(endogenous_vars = excluded_instruments)
```

- **endogenous_vars**: Variables treated as endogenous (instrumented)
- **excluded_instruments**: Instruments not included in the structural equation
- **exogvars**: Exogenous regressors (included in both stages)

The full instrument set consists of the excluded instruments plus all exogenous regressors.

## Algorithm

`civreghdfe` implements 2SLS estimation in four steps:

### Step 1: HDFE Absorption
Absorb fixed effects from all variables (y, endogenous X, exogenous X, instruments Z) using the conjugate gradient (CG) solver with symmetric Kaczmarz iteration.

### Step 2: First Stage
Regress each endogenous variable on all instruments (excluded + exogenous) using demeaned data.

### Step 3: Second Stage (2SLS)
Compute the 2SLS estimator on demeaned data:
```
β = (X'Pz X)⁻¹ X'Pz y
```
where Pz = Z(Z'Z)⁻¹Z' is the projection onto the instrument space.

### Step 4: Variance Estimation
Compute residuals using original X (not projected X) for proper inference. Apply the appropriate sandwich estimator.

## First-Stage Statistics

With the `first` option, `civreghdfe` reports:
- First-stage F-statistics for each endogenous variable
- Tests instrument strength (rule of thumb: F > 10)

Computed from the partial R-squared of excluded instruments in first-stage regressions.

## Variance-Covariance Estimation

| VCE Type | Description |
|----------|-------------|
| `unadjusted` | Standard VCE assuming homoskedasticity |
| `robust` | Heteroskedasticity-robust (HC1) |
| `cluster clustvar` | Cluster-robust, clustered on `clustvar` |

All VCE computations use the proper 2SLS sandwich formula with original endogenous regressors (not first-stage fitted values).

## Stored Results

`civreghdfe` stores the following in `e()`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(df_r)` | Residual degrees of freedom |
| `e(df_a)` | Absorbed degrees of freedom |
| `e(K)` | Number of regressors |
| `e(K_endog)` | Number of endogenous regressors |
| `e(K_exog)` | Number of exogenous regressors |
| `e(K_iv)` | Number of instruments (including exogenous) |
| `e(G)` | Number of absorbed FE groups |
| `e(N_clust)` | Number of clusters (if clustered) |
| `e(F_first1)` | First-stage F-stat for 1st endogenous var |
| `e(F_first2)` | First-stage F-stat for 2nd endogenous var |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `civreghdfe` |
| `e(cmdline)` | Command as typed |
| `e(depvar)` | Dependent variable |
| `e(endogvars)` | Endogenous variables |
| `e(exogvars)` | Exogenous variables |
| `e(instruments)` | Excluded instruments |
| `e(absorb)` | Absorbed fixed effects |
| `e(vcetype)` | VCE type |
| `e(clustvar)` | Cluster variable (if clustered) |

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix |

## Identification Requirements

For valid 2SLS estimation:
1. **Order condition**: Number of excluded instruments ≥ number of endogenous variables
2. **Rank condition**: Instruments must be correlated with endogenous variables
3. **Exclusion restriction**: Instruments must be uncorrelated with structural error

## Technical Notes

### HDFE Absorption
- Uses same CG algorithm as `creghdfe`
- Iteratively projects out factor means
- Convergence controlled by `tolerance()` and `maxiter()`

### Degrees of Freedom
- Absorbed fixed effects consume degrees of freedom
- Proper adjustment for clustered standard errors

### Overidentification
- When instruments > endogenous variables, the model is overidentified
- Overidentification tests not yet implemented

## Performance

`civreghdfe` achieves its speed through the following tricks:

**HDFE Absorption (shared with `creghdfe`):**
- **CG solver with symmetric Kaczmarz transform**: Forward + backward sweep per iteration doubles convergence rate
- **Pre-computed inverse group counts**: Multiplies instead of divides in the inner demeaning loop
- **Counting sort for CSR construction**: O(N+L) sparse index builds
- **Thread-local buffers**: Eliminates contention during parallel demeaning

**2SLS Estimation:**
- **Cholesky-based normal equations**: Exploits symmetric positive-definiteness of Z'Z and X'P_Z X for efficient inversion
- **Fused project-subtract**: Single-pass FE projection reduces memory traffic
- **K-way unrolled dot product**: 8/16-element unrolling for inner products in both stages

**Data I/O and Infrastructure:**
- **SIMD bulk load/store**: AVX2/SSE2/NEON with 16x unrolling via `ctools_data_io`
- **OpenMP parallelization**: Parallel matrix operations and data loading
- **Cache-line aligned allocations**: 64-byte boundaries prevent false sharing
- **Persistent thread pool**: Threads reused across first stage, second stage, and HDFE iterations

Typical speedup: 5-10x over `ivreghdfe`.

## Comparison with ivreghdfe

| Feature | `ivreghdfe` | `civreghdfe` |
|---------|-------------|--------------|
| Implementation | Stata + Mata | C with OpenMP |
| HDFE Algorithm | Alternating projections | CG with Kaczmarz |
| Parallelization | Limited | Yes |
| Typical Speedup | 1x (baseline) | 5-10x |

## See Also

- [ctools Overview](../README.md)
- [creghdfe](README_creghdfe.md) - OLS regression with HDFE
- ivreghdfe - Original implementation (if installed)
- [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html) - Alternative IV estimation

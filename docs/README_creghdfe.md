# creghdfe

High-performance C-accelerated high-dimensional fixed effects regression.

## Overview

`creghdfe` is a drop-in replacement for `reghdfe` that performs linear regression with multiple high-dimensional fixed effects. All computation is performed in optimized C code using a conjugate gradient (CG) solver with Kaczmarz-style sweeping, achieving up to 10x speedup over pure Stata/Mata implementations.

## Syntax

```stata
creghdfe depvar indepvars [if] [in] [weight], absorb(varlist) [options]
```

## Options

### Model Options
| Option | Description |
|--------|-------------|
| `absorb(varlist)` | **Required.** Categorical variables representing fixed effects to absorb |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce(vcetype)` | Variance-covariance estimator |

Where `vcetype` can be:
- `unadjusted` or `ols` - Standard VCE assuming homoskedasticity
- `robust` - Heteroskedasticity-robust (HC1) standard errors
- `cluster varlist` - Cluster-robust standard errors

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information and timing |
| `timeit` | Display total elapsed time |

### Weights
| Weight Type | Description |
|-------------|-------------|
| `aweight` | Analytic weights |
| `fweight` | Frequency weights |
| `pweight` | Probability weights (forces robust VCE) |

## Examples

```stata
* Basic regression with firm and year fixed effects
creghdfe y x1 x2, absorb(firm year)

* With clustered standard errors
creghdfe y x1 x2, absorb(firm year) vce(cluster firm)

* With robust standard errors
creghdfe y x1 x2, absorb(firm year) vce(robust)

* With analytic weights
creghdfe y x1 x2 [aw=weight], absorb(firm year)

* With frequency weights
creghdfe y x1 x2 [fw=freq], absorb(firm year)

* With probability weights (automatically uses robust VCE)
creghdfe y x1 x2 [pw=sampwt], absorb(firm year)

* Verbose output showing timing breakdown
creghdfe y x1 x2, absorb(firm year) verbose

* Multiple absorbed dimensions
creghdfe y x1 x2, absorb(firm year industry)
```

## Performance

`creghdfe` achieves its speed through:

### Conjugate Gradient Solver
- Uses CG iteration with Kaczmarz-style symmetric sweeping
- Converges in O(sqrt(condition number)) iterations
- Much faster than direct projection for large fixed effects

### Parallel Computation
- OpenMP parallelization for matrix operations
- Parallel data loading with 8-way loop unrolling
- Efficient memory access patterns

### Iterative Demeaning
- Fixed effects are absorbed via iterative mean removal
- Convergence is fast for typical panel data structures

## Algorithm

The fixed effects absorption works as follows:

1. **Data Loading**: Load dependent and independent variables into C memory
2. **Demeaning**: Iteratively demean all variables within each absorbed factor
3. **OLS Estimation**: Solve the normal equations on demeaned data
4. **VCE Computation**: Calculate standard errors with appropriate adjustments

The CG solver minimizes:
```
||y - Xβ - Σ αᵢ||²
```
where αᵢ represents the fixed effects for each absorbed variable.

## Stored Results

`creghdfe` stores results in `e()` similar to `reghdfe`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(df_r)` | Residual degrees of freedom |
| `e(df_m)` | Model degrees of freedom |
| `e(df_a)` | Degrees of freedom absorbed by FEs |
| `e(r2)` | R-squared |
| `e(r2_a)` | Adjusted R-squared |
| `e(r2_within)` | Within R-squared |
| `e(rss)` | Residual sum of squares |
| `e(rmse)` | Root mean squared error |
| `e(N_clust)` | Number of clusters (if clustered) |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `creghdfe` |
| `e(depvar)` | Dependent variable name |
| `e(absorb)` | Absorbed fixed effects |
| `e(vce)` | VCE type |
| `e(clustvar)` | Cluster variable (if clustered) |

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix |

## Technical Notes

- The absorbed fixed effects consume degrees of freedom
- Singleton groups (groups with only one observation) are dropped
- Multicollinear fixed effects are automatically handled
- The constant term is absorbed with the fixed effects

## Comparison with reghdfe

| Feature | `reghdfe` | `creghdfe` |
|---------|-----------|------------|
| Implementation | Stata + Mata | C with OpenMP |
| Algorithm | Method of Alternating Projections | CG with Kaczmarz |
| Parallelization | Limited | Yes |
| Typical Speedup | 1x (baseline) | 5-10x |
| Memory | Higher | Optimized |

## Degrees of Freedom

The residual degrees of freedom is computed as:
```
df_r = N - K - df_absorbed + redundant_FEs
```
where redundant FEs are fixed effects that are collinear with others.

## See Also

- [ctools Overview](../README.md)
- [reghdfe](http://scorreia.com/software/reghdfe/) - Original Stata implementation
- [cqreg](README_cqreg.md) - Quantile regression with HDFE
- [civreghdfe](README_civreghdfe.md) - IV regression with HDFE

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

### Speedup Tricks

**CG Solver and HDFE Absorption:**
- **Symmetric Kaczmarz transform**: Forward + backward sweep per CG iteration doubles convergence rate compared to a one-directional sweep
- **Pre-computed inverse group counts**: Multiplies by `1.0/count` instead of dividing, since division is 5-20x slower than multiplication on modern CPUs
- **Counting sort for CSR index construction**: Builds compressed sparse row indices in O(N+L) instead of O(N log N), where L is the number of FE levels
- **Two-pass singleton detection**: First pass counts group sizes, second pass marks singletons—avoids expensive repeated scans
- **Iterative singleton removal**: Handles connected components where dropping one singleton creates new singletons

**Linear Algebra:**
- **K-way unrolled dot product**: `ctools_dot_unrolled` processes 8 or 16 elements per loop iteration, improving instruction-level parallelism and reducing loop overhead
- **Fused project-subtract**: Combines FE projection and subtraction in a single pass over the data, halving memory traffic
- **Cholesky decomposition for (X'X)^-1**: Exploits symmetric positive-definiteness for efficient normal equation solving
- **OpenMP SIMD pragmas**: Compiler hints trigger auto-vectorization of inner loops

**Data I/O:**
- **SIMD bulk load/store**: AVX2/SSE2 (x86) and NEON (ARM64) for parallel variable loading with 16x unrolling
- **Thread-local buffers**: Each OpenMP thread has its own accumulator, eliminating contention during parallel demeaning
- **Cache-line aligned allocations**: 64-byte boundaries prevent false sharing between threads
- **Persistent thread pool**: Threads are reused across CG iterations

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

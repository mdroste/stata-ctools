# cbinscatter

High-performance C-accelerated binned scatter plots.

## Overview

`cbinscatter` is a high-performance replacement for `binscatter` that uses a C plugin for all data computations. It creates binned scatter plots by dividing the x variable into equal-sized quantile bins and plotting the mean of y within each bin. All heavy computation (data loading, residualization, bin computation, and line fitting) is performed in optimized C code with OpenMP parallelization.

## Syntax

```stata
cbinscatter yvar xvar [if] [in] [weight] [, options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `nquantiles(#)` | Number of bins (default: 20, range: 2-1000) |
| `controls(varlist)` | Control variables to partial out via OLS |
| `absorb(varlist)` | Fixed effects to absorb via HDFE |
| `by(varname)` | Create separate series by group |

### Line Fitting Options
| Option | Description |
|--------|-------------|
| `linetype(type)` | Fit line type: `none`, `linear` (default), `qfit`/`quadratic`, `cubic` |

### Data Options
| Option | Description |
|--------|-------------|
| `discrete` | Treat x as discrete (one bin per unique value) |
| `genxq(varname)` | Generate bin assignment variable (not yet implemented) |
| `savedata(filename)` | Save bin data to Stata file |

### Graph Options
| Option | Description |
|--------|-------------|
| `nograph` | Suppress graph output |
| `title(string)` | Graph title |
| `ytitle(string)` | Y-axis title (default: variable name) |
| `xtitle(string)` | X-axis title (default: variable name) |
| `legend(string)` | Legend options |
| `colors(string)` | Colors for series |
| *twoway_options* | Additional options passed to `twoway` |

### Reporting Options
| Option | Description |
|--------|-------------|
| `reportreg` | Report underlying regression |
| `verbose` | Display progress information |
| `timeit` | Display timing breakdown |

### Weights
Supports `aweight`, `fweight`, `pweight`, and `iweight`.

## Examples

```stata
* Basic binned scatter plot
sysuse auto
cbinscatter price mpg

* With 10 bins and no fit line
cbinscatter price mpg, nquantiles(10) linetype(none)

* Control for other variables
cbinscatter price mpg, controls(weight length)

* With fixed effects
webuse nlswork, clear
cbinscatter ln_wage tenure, absorb(idcode year)

* Separate series by group
sysuse auto
cbinscatter price mpg, by(foreign) nquantiles(10)

* Quadratic fit with custom titles
cbinscatter price mpg, linetype(qfit) title("Price vs Fuel Efficiency") ///
    ytitle("Price (USD)") xtitle("Miles per Gallon")

* With weights
cbinscatter price mpg [aw=weight], nquantiles(15)

* Save bin data without displaying graph
cbinscatter price mpg, nograph savedata(mybins)

* Cubic fit with verbose timing
cbinscatter price mpg, linetype(cubic) verbose timeit
```

## Performance

`cbinscatter` achieves dramatic speedups through:

### Histogram-Based Binning
- Uses O(N) histogram algorithm instead of O(N log N) sorting
- Sequential memory access for cache efficiency
- 40x+ speedups for bin computation alone

### Parallel Data Loading
- 8-way loop unrolling for efficient data transfer
- OpenMP parallelization where applicable

### Single-Pass Statistics
- Bin means computed in single pass
- Linear and quadratic fits use closed-form solutions

### Optimized Residualization
- Controls partialled out via efficient OLS
- Fixed effects absorbed using `creghdfe`'s HDFE algorithm

### Performance Benchmark
On a dataset with 25 million observations, `cbinscatter` completes in under 1 second including graph generation.

## Timing Breakdown

With `timeit`, you'll see time spent in:
- **Load**: Loading data from Stata
- **Residualize**: Partialling out controls/fixed effects
- **Bin computation**: Computing bin means
- **Line fitting**: Fitting trend line
- **Graph**: Generating the plot

## Stored Results

`cbinscatter` stores the following in `e()`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations used |
| `e(N_dropped)` | Observations dropped (missing values) |
| `e(nquantiles)` | Number of bins |
| `e(num_groups)` | Number of by-groups |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `cbinscatter` |
| `e(cmdline)` | Command as typed |
| `e(depvar)` | Y variable |
| `e(xvar)` | X variable |
| `e(controls)` | Control variables |
| `e(absorb)` | Absorbed fixed effects |
| `e(by)` | By variable |
| `e(linetype)` | Fit line type |

### Matrices
| Result | Description |
|--------|-------------|
| `e(bindata)` | Bin statistics: by_group, bin_id, x_mean, y_mean, x_se, y_se, n_obs |
| `e(coefs)` | Fit line coefficients (if linetype specified) |
| `e(fit_stats)` | R-squared and other fit statistics |

## How Binning Works

1. **Residualization** (if controls or absorb specified):
   - OLS residualization for controls
   - HDFE absorption for fixed effects

2. **Bin Assignment**:
   - Observations sorted into equal-sized quantile bins
   - Uses histogram approach for O(N) complexity

3. **Bin Statistics**:
   - Mean of x and y computed within each bin
   - Standard errors optionally computed

4. **Line Fitting**:
   - Fit computed from microdata (not bin means)
   - Supports linear, quadratic, and cubic fits

## Technical Notes

- Missing values in x, y, controls, or absorb variables drop observations
- With `by()`, each group gets independent bins and fit lines
- The `discrete` option creates one bin per unique x value
- Fit line coefficients are from the residualized data

## Comparison with binscatter

| Feature | `binscatter` | `cbinscatter` |
|---------|--------------|---------------|
| Implementation | Stata + Mata | C with OpenMP |
| Binning Algorithm | Sort-based | Histogram-based |
| Binning Complexity | O(N log N) | O(N) |
| Typical Speedup | 1x (baseline) | 10-100x |
| HDFE Support | Via reghdfe | Native |

## See Also

- [ctools Overview](../README.md)
- [creghdfe](README_creghdfe.md) - Regression with HDFE (used internally)
- binscatter - Original Stata implementation (if installed)

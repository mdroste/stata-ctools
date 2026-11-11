# cbinscatter

Create binned scatter plots.

## Overview

`cbinscatter` groups x values into quantile bins and plots the mean of y against the mean of x in each bin. C code performs the calculations; Stata draws the graph.

With `method(binsreg)`, observations in singleton fixed-effect groups remain in
the bins and in the sample used to center adjusted means. Comparisons with
`binsreg` should use `reghdfeopt(keepsingletons)` and the same bin partition;
dropping singletons can shift the adjusted means even when bin contrasts agree.
The validation suite supplies explicit knots between `xtile` groups and saves
reference dots in double precision, avoiding the rounding of cutpoints and saved
coordinates in the reference package. It retains the standard seven-significant-figure
comparison threshold.

## Syntax

```stata
cbinscatter yvar [yvar ...] xvar [if] [in] [weight] [, options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `nquantiles(#)` (`n()`) | Number of bins (default: 20, at least 2; all by() groups' bins must fit in `c(max_matdim)` rows). binscatter stops with r(198) above 100 quantile bins; cbinscatter forms them |
| `controls(varlist)` | Control variables to partial out via OLS |
| `absorb(varlist)` | Fixed effects to absorb via HDFE (numeric or string; any codes) |
| `by(varname)` | Create separate series by group |
| `method(string)` | `classic` (default) residualizes x and y; `binsreg` bins raw x and adjusts bin means jointly with controls |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |

### Line Fitting Options
| Option | Description |
|--------|-------------|
| `linetype(type)` | Fit line type: `none`, `linear` (default), `qfit`/`quadratic`, `cubic` |

### Data Options
| Option | Description |
|--------|-------------|
| `discrete` | Treat x as discrete (one bin per unique value) |
| `xq(varname)` | Bins given by the values (positive integers) of a variable, as binscatter's `xq()` |
| `genxq(newvar)` | Store each observation's quantile bin, numbered as binscatter's `genxq()` numbers it |
| `medians` | Plot bin medians (`_pctile`'s, weighted) instead of means |
| `noaddmean` | Keep the residualized x and y without their means (with `controls()`/`absorb()`) |
| `rd(numlist)` | Fit lines separately on each side of these x values; dashed lines at them |
| `randvar(varname)`, `randcut(#)`, `randn(#)` | Quantile cutpoints from a random subsample, as binscatter's `fastxtile` (the same draws under the same seed) |
| `savedata(filename)` | Without an extension: binscatter's `filename.csv` and `filename.do` (classic method), plus the bins in `filename.dta`; with an extension, or with `method(binsreg)`, the `.dta` file alone |
| `replace` | Let `savedata()` and `savegraph()` overwrite files |

`medians`, `noaddmean`, `rd()`, `xq()` and the random-subsample options are binscatter's (classic method) and
are refused with `method(binsreg)`. All of them work with several y variables, `by()`, `controls()` and `absorb()`.

### Graph Options
| Option | Description |
|--------|-------------|
| `nograph` | Suppress graph output |
| `title(string)` | Graph title |
| `ytitle(string)` | Y-axis title (default: variable name) |
| `xtitle(string)` | X-axis title (default: variable name) |
| `legend(string)` | Legend options |
| `colors(string)` | Colors for series |
| `mcolors(string)`, `lcolors(string)` | Marker and line colors for series, as in binscatter |
| `savegraph(filename)` | `graph save` (.gph or no extension) or `graph export` the graph; needs `replace` to overwrite |
| `msymbols(string)` (`m()`) | Marker symbols for series |
| `mlabels(string)` | Text labels for the plotted series |
| *twoway_options* | Additional options passed to `twoway` |

### binscatter legacy options
`nbins(#)`, `create_xq`, `x_q(varname)`, `symbols()`, `unique()` and `linetype(noline)` work as in
binscatter 7.02, with its notes. With both an option and its legacy form binscatter shows its message
and stops with r(0); cbinscatter shows the message and exits with r(198).

### Differences from binsreg
`method(binsreg)` computes binsreg's binned means on fixed quantile bins of the raw x. binsreg was not
installed to verify against, so these binsreg features are not provided: IMSE-optimal bins
(`nbins()` selection, `binsmethod()`, `nbinsrot()`, `binspos()`, `pselect()`/`sselect()`), `dots(p s)`
with p > 0, `line()`, `ci()`, `cb()`, `deriv()`, `at()`, `vce()`, `level()`, `nsims()`; `polyreg(p)` above 3
(at most `linetype(cubic)`); `samebinsby`, `bycolors()`, `bysymbols()`, `bylpatterns()` (each group is
binned on its own; series take `colors()`, `mcolors()`, `lcolors()`, `msymbols()`); `masspoints()`,
`dfcheck()`, `noplot` (use `nograph`), `plotxrange()`/`plotyrange()` (use `xscale()`/`yscale()`),
`reghdfeopt()`, binsreg's `savedata()` layout, and binsreg's `e()` results. binsreg's positional controls
(`binsreg y x w1 w2`) go in `controls(w1 w2)`: more than two variables are several y and one x, as in
binscatter.

### Reporting Options
| Option | Description |
|--------|-------------|
| `reportreg` | Report underlying regression |
| `verbose` | Display progress information |

### Weights
Supports `aweight`, `fweight`, `pweight`, and `iweight`.

## Examples

```stata
* Basic binned scatter plot
sysuse auto, clear
cbinscatter price mpg

* With 10 bins and no fit line
cbinscatter price mpg, nquantiles(10) linetype(none)

* Control for other variables
cbinscatter price mpg, controls(weight length)

* With fixed effects
webuse nlswork, clear
cbinscatter ln_wage tenure, absorb(idcode year)

* Separate series by group
sysuse auto, clear
cbinscatter price mpg, by(foreign) nquantiles(10)

* Quadratic fit with custom titles
cbinscatter price mpg, linetype(qfit) title("Price vs Fuel Efficiency") ///
    ytitle("Price (USD)") xtitle("Miles per Gallon")

* With weights
cbinscatter price mpg [aw=weight], nquantiles(15)

* Save bin data without displaying graph
cbinscatter price mpg, nograph savedata(mybins)

* Cubic fit with verbose timing
cbinscatter price mpg, linetype(cubic) verbose
```

## Performance

`cbinscatter` uses the following methods to compute bins and summary statistics:

### Implementation details

**Binning Algorithm:**
- **Exact quantile binning in O(N)**: Sorts the x values with a radix sort and uses the weighted percentile definitions of `xtile` and `_pctile` to determine bin boundaries
- **Direct bin assignment from sorted order**: When data is sorted, bin boundaries are read off directly without binary search
- **O(N) weighted quantile computation**: Accumulates cumulative weights in a single pass to find quantile boundaries

**Statistics Computation:**
- **Single-pass accumulation**: Mean, variance, and standard error for each bin are computed in one pass over the data, avoiding repeated traversals
- **Contiguous allocation block**: All per-bin accumulators (x_sum, y_sum, x2_sum, count, etc.) are allocated in a single cache-friendly memory block
- **Fast path for unweighted no-SE case**: A specialized loop skips weight and SE computation when neither is needed, reducing per-observation work

**Residualization:**
- **Iterative projections**: `cbinscatter_resid.c` absorbs fixed effects by repeated projection sweeps; it does not use the shared conjugate-gradient solver in `creghdfe`.
- **Fused project-subtract**: Combines projection and subtraction in a single memory pass

**Data I/O:**
- **SIMD bulk load**: AVX2/SSE2/NEON accelerated with 16x loop unrolling via `ctools_data_io`
- **Overflow-safe allocation**: `ctools_safe_mul_size` checks prevent silent overflow on datasets with many bins or groups

### Performance Benchmark
Use the validation benchmarks to measure a particular dataset and configuration. Record the source revision, dimensions, options, hardware, and reference package versions. Runtime depends on the dataset, options, hardware, and thread count.

## Timing Breakdown

With `verbose`, you'll see time spent in:
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
| `e(discrete)` | 1 if the bins are the values of x |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `cbinscatter` |
| `e(cmdline)` | Command as typed |
| `e(depvar)` | Y variable (variables) |
| `e(xvar)` | X variable |
| `e(controls)` | Control variables |
| `e(absorb)` | Absorbed fixed effects |
| `e(by)` | By variable |
| `e(linetype)` | Fit line type |
| `e(rd)` | With `rd()`, the break points |
| `e(graphcmd)` | The `twoway` command that drew the graph |

### Matrices
| Result | Description |
|--------|-------------|
| `e(bindata)` | Bin statistics: by_group, bin_id, x_mean, y_mean, x_se, y_se, n_obs |
| `e(coefs)` | Fit line coefficients per by-group (and `rd()` interval): constant, x, x^2, x^3 (if linetype specified) |
| `e(fit_stats)` | R-squared and number of observations of each fit |
| `e(rdintervals)` | With `rd()`, each interval's bounds (`gt`, `lt_eq`), as in binscatter |
| `e(y#_coefs)` | binscatter's fit coefficients of the #th y: x terms from the highest power, `_cons`; rows `by#:rd#` |
| `e(byvalues)` | With a numeric `by()`, its values (rows `by1`, `by2`, ...), as in binscatter |

### Functions
| Result | Description |
|--------|-------------|
| `e(sample)` | The estimation sample, as in binscatter |

## How Binning Works

1. **Residualization** (if controls or absorb specified):
   - OLS residualization for controls
   - HDFE absorption for fixed effects

2. **Bin Assignment**:
   - Observations sorted into equal-sized quantile bins (tied x values share a bin)
   - Radix sort of x for O(N) complexity; cutpoints match `xtile`

3. **Bin Statistics**:
   - Mean of x and y computed within each bin
   - Standard errors optionally computed

4. **Line Fitting**:
   - Fit computed from microdata (not bin means)
   - Supports linear, quadratic, and cubic fits
   - Powers of centered, scaled x; redundant terms are omitted, as `regress` does

## Technical Notes

- Missing values in x, y, controls, or absorb variables drop observations
- `absorb()` values only label groups: negative, fractional, large (e.g. 11-digit) and string codes are allowed
- With `by()`, groups share the pooled-sample quantile bins (as `binscatter` does; `method(binsreg)` bins each group separately) and get their own fit lines
- The `discrete` option creates one bin per unique x value; with `method(binsreg)` those means are adjusted for controls/absorb like quantile bins
- `method(classic)` with controls/absorb fits the line to the residualized data with the sample means added back (binscatter's convention), so `e(coefs)` matches binscatter's `e(y1_coefs)` (which lists `_cons` last)
- `method(binsreg)` fits the polynomial jointly with the controls and absorbed effects and evaluates it at the controls' means (binsreg's `polyreg()` convention)

## Comparison with binscatter

| Feature | `binscatter` | `cbinscatter` |
|---------|--------------|---------------|
| Implementation | Stata + Mata | C with OpenMP |
| Binning Algorithm | Sort-based | Radix sort, same cutpoints |
| Binning Complexity | O(N log N) | O(N) |
| HDFE Support | Via reghdfe | Native |

## See Also

- [ctools Overview](../../README.md)
- [creghdfe](creghdfe.md) - Regression with HDFE
- binscatter - Original Stata implementation (if installed)

## Controls and weights

Weight expressions such as `[aw=1+abs(z)]` are evaluated once on the candidate
sample. Missing weights are excluded; selected weights must be positive, and
frequency weights must be integers. `controls()` accepts factor terms and
interactions, which the command expands into indicator and interaction variables
before estimation. For
example, `controls(i.category c.z#i.category)` is supported. The two plotted
variables must be numeric variables.

`method(binsreg)` omits redundant control columns when solving the adjustment,
including with absorbed effects. Adding a duplicate control leaves fitted bin
means unchanged. Failure to solve the adjustment returns an error.

# Compatibility

Each ctools command aims to work like the command it replaces, and for everyday use the syntax is the same. Not every option is implemented, some commands add options of their own, and the estimators can differ from the originals in solver details. This page covers those differences. For full syntax, see `help <command>` in Stata.

In the Weights column, `aw`, `fw`, `pw`, and `iw` are Stata weight types. A dash means the command doesn't accept weights.

| Command | Replaces | Weights | Notes |
|---|---|---|---|
| `cimport` | `import` | — | FRED, Haver, and StataNow's Parquet import aren't supported. See [format-specific differences](IO_PARITY.md). |
| `cexport` | `export` | — | See [format-specific differences](IO_PARITY.md). |
| `cdistinct` | `distinct` | — | Uses the original `in` range for each variable; distinct 1.2.1 can reselect rows after sorting. strL text is limited to 2,045 bytes. |
| `chdfe` | `hdfe` | aw/fw/pw | Categorical interactions and heterogeneous slopes are supported. Solver iteration paths may differ. Conjugate gradient requires a symmetric transform. See [details](commands/chdfe.md). |
| `csort` | `sort` | — | Stable by default. Rejects datasets containing any strL variables, even when already sorted. Adds `stream(#)` and `nosortedby` options. |
| `cmerge` | `merge` | — | Key variables must be the same type (string or numeric) in both datasets, even with `force`. |
| `csample` | `sample` | — | Uses `csample, count(20)` instead of `sample 20, count`. Samples need not match `sample` with the same seed. |
| `cbsample` | `bsample` | — | — |
| `cencode` | `encode` | — | Accepts several variables at once and adds `replace`. |
| `cdecode` | `decode` | — | Accepts several variables at once and adds `replace`. |
| `cdestring` | `destring` | — | strL values are limited to 2,045 bytes. |
| `cinterpolate` | `ipolate`, `cipolate` | — | Linear by default; `cubic` selects local four-point cubic interpolation without extrapolation. Adds `by()` and keeps row order. Cubic anchors use distinct x values after averaging ties; all extended missings are excluded. See [usage](USAGE.md#interpolation-cinterpolate) for migration from the former ctools name. |
| `csplit` | `split` | — | — |
| `crangejoin` | `rangejoin` 1.1.3 | — | Numeric and str# variables only. Output is limited to 2,147,483,647 observations. |
| `cwinsor` | `winsor2`, `gstats winsor` | — | — |
| `cbinscatter` | `binscatter`, `binsreg` | aw/fw/pw/iw | `method(binsreg)` omits confidence bands and tests. It keeps singleton fixed-effect groups; comparisons with `binsreg` require the same bins and `reghdfeopt(keepsingletons)` there. |
| `crangestat` | `rangestat` | — | User-written Mata functions aren't supported. |
| `creghdfe` | `reghdfe` | aw/fw/pw | `absorb()` takes variables, `i.var`, and `var1#var2`, but not heterogeneous slopes (`c.var`), `##`, or `newvar=` terms. Clustering is one-way only. Not every `reghdfe` solver or degrees-of-freedom option is available. |
| `civreghdfe` | `ivreghdfe` | aw/fw/pw | `absorb()` doesn't take heterogeneous slopes, and `center` isn't supported. Some estimator and VCE combinations aren't available; see the help file. |
| `cqreg` | `qreg` | — | Adds clustered standard errors. Different solver and density estimator can produce slightly different coefficients and standard errors. Errors on nonconvergence. |
| `cpsmatch` | `psmatch2` | — | ATT for one outcome only; no ATE/ATU or `ai()`. Matching without replacement requires `neighbor(1)`. |
| `cppmlhdfe` | `ppmlhdfe` | aw/fw/pw | Undocumented debugging options aren't available. The iteration path can differ even when final estimates agree. See [cppmlhdfe details](commands/cppmlhdfe.md). |

## cpsmatch standard errors

`cpsmatch` reports the same ATT standard error as `psmatch2` does by default:

```
SE = sqrt(s1²/n1 + s0² · Σ w0² / n1²)
```

Here `n1` is the number of matched treated observations, `s1²` and `s0²` are the unweighted outcome variances among matched treated observations and the controls they're matched to, and `w0` are the controls' matching weights. Squaring the weights accounts for controls used more than once and for fractional matches. Like `psmatch2`'s default, this ignores the uncertainty from estimating the propensity score.

The Unmatched row also follows `psmatch2`: its standard error comes from an OLS regression of the outcome on treatment.

## String variables

The shared C transfer layer handles numeric variables and fixed-width strings up to `str2045`. strL support depends on the command: `cencode`, `cdecode`, `csplit`, and `cmerge` use native Stata paths where needed, while `csort` and `crangejoin` reject strL. `cdestring` limits strL text to 2,045 bytes. `cimport` and `cexport` have dedicated long-text paths for delimited, Excel, and statistical formats, and handle shapefile binary headers. See the command help for format-specific limits.

## When a command fails

`cmerge`, `crangejoin`, `cwinsor`, `crangestat`, `creghdfe`, and `cimport` undo their changes if they fail or you press Break, without writing a copy of the dataset to disk. Two options change this:

- `snapshot` runs `preserve` on the whole dataset first. It's slower on large data, since Stata writes the snapshot to disk once it exceeds `c(max_preservemem)`.
- `fast` keeps no copy at all. It saves time and memory, but a failure or Break can leave the data partly changed.

On Stata 14 and 15, which don't have frames, `cimport`, `crangejoin`, and some `cmerge` merges use `preserve` instead.

Like Stata's own commands, `cimport` clears the data before reading SAS, SPSS, XPORT, dBase, and shapefiles, so a failed import of those formats leaves an empty dataset unless you specify `snapshot`.

If any command says the data in memory may be incomplete or partly modified, reload your data before you continue.

New variables created by `creghdfe` (residuals, group IDs, saved fixed effects), `crangestat`, and `cencode` must have names that aren't already in use.

## Estimation errors

- **One cluster.** A clustered covariance matrix isn't defined with a single cluster. `creghdfe` and `civreghdfe` exit with error 498. `cppmlhdfe` exits with error 459, since it needs at least two clusters in every cluster dimension.
- **No convergence.** If the fixed-effect projection fails in `creghdfe`, `civreghdfe`, or `cppmlhdfe`, or `cppmlhdfe` runs out of iterations, the command exits with error 430 instead of reporting estimates.

## Performance

Speedups depend on the number of observations and variables, variable types, options, hardware, and thread count. In some cases ctools is slower than the command it replaces. See the [benchmarks](performance/README.md) for measured results.

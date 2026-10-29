# Compatibility

Each ctools command aims to work like the command it replaces, and for everyday use the syntax is the same. Not every option is implemented, some commands add options of their own, and the estimators can differ from the originals in solver details. This page covers those differences. For full syntax, see `help <command>` in Stata.

In the Weights column, `aw`, `fw`, `pw`, and `iw` are Stata weight types. A dash means the command doesn't accept weights.

| Command | Replaces | Weights | Notes |
|---|---|---|---|
| `cimport` | `import` | — | Reads delimited text, XLS/XLSX, dBase, SAS (including catalogs), SAS XPORT v5/v8, SPSS SAV/ZSAV, and shapefiles. Handles Unicode and code-page encodings, custom delimiters, and locale-specific number formats. Long text and shapefile binary headers are stored as strL. FRED, Haver, and StataNow's Parquet import aren't supported. See [cimport details](README_cimport.md). |
| `cexport` | `export` | — | Writes delimited text, XLS/XLSX, dBase, SAS, SAS XPORT v5/v8, SPSS SAV, and shapefiles. Can update sheets in existing XLS and XLSX workbooks. Output goes to a temporary file first and replaces the destination only once it's complete, so a failed export doesn't leave a half-written file. See [cexport details](README_cexport.md). |
| `csort` | `sort` | — | Stable ascending sort. Doesn't take `if` or `in`, and can't sort on strL variables. `stream(#)` reduces memory use on wide datasets. `nosortedby` skips marking the data as sorted, which saves time when nothing afterward relies on Stata knowing the sort order. |
| `cmerge` | `merge` | — | Supports the same merge types and main options. Checks that the key is unique (for 1:1, m:1, and 1:m) before changing anything, and restores the master data if the merge fails. Key variables must be the same type (string or numeric) in both datasets, even with `force`. |
| `csample` | `sample` | — | Samples without replacement by percentage (`csample 10`) or count (`csample, count(20)`, where `sample` uses `sample 20, count`). Supports `by()`. Results are reproducible with `set seed`, but they won't match the observations `sample` picks with the same seed. |
| `cbsample` | `bsample` | — | Bootstrap resampling of observations or clusters, optionally within strata. The number you give is the number of draws per stratum, counted in observations or clusters. `weight(varname)` fills an existing numeric variable with frequency weights instead of resampling. |
| `cencode` | `encode` | — | Encodes several variables at once and adds a `replace` option. Works with str# and strL variables. In a few edge cases, such as existing value labels that contain extended missing values or very long text, it hands the work to `encode`. |
| `cdecode` | `decode` | — | Gives the same values and storage types as `decode`. Decodes several variables at once and adds a `replace` option that converts all of them or none. Labels longer than 2,045 bytes, and `maxlength()` cuts through multibyte characters, are handled by `decode` itself, so `threads()` doesn't speed those up. |
| `cdestring` | `destring` | — | Converts str# variables and strL values up to 2,045 bytes. Supports `ignore()`, `force`, and the usual output options. |
| `cipolate` | `ipolate` | — | Linear interpolation and extrapolation (`epolate`), with `by()` or the `by:` prefix and `if`/`in`. Keeps the original row order. Not related to the cubic-interpolation command of the same name on SSC. |
| `csplit` | `split` | — | Supports `parse()`, `notrim`, `limit()`, `if`/`in`, and `split`'s `destring` options. strL variables are split by Stata's own `split`. New variables are added only if the whole command succeeds. |
| `crangejoin` | `rangejoin` 1.1.3 | — | Joins on inclusive intervals, with `by()`, `keepusing()`, `prefix()`, `suffix()`, and `all`. Keeps unmatched master observations, and matches stay in the using data's order. Numeric and str# variables only. Output is limited to 2,147,483,647 observations. |
| `cwinsor` | `winsor2`, `gstats winsor` | — | Winsorizes or trims, optionally within groups, into new or replaced variables. New variables are missing outside `if`/`in`; with `replace`, those observations keep their original values. |
| `cbinscatter` | `binscatter`, `binsreg` | aw/fw/pw/iw | Plotted variables must be numeric; controls can be factor variables. Redundant controls are dropped. `method(binsreg)` follows `binsreg`'s approach to binning and adjustment but doesn't compute its confidence bands or tests. It keeps singleton fixed-effect groups, so to compare with `binsreg`, use the same bins and `reghdfeopt(keepsingletons)` there. |
| `crangestat` | `rangestat` | — | Built-in statistics over intervals, with `by()` and `excludeself`. User-written Mata functions aren't supported. |
| `creghdfe` | `reghdfe` | aw/fw/pw | `absorb()` takes variables, `i.var`, and `var1#var2`, but not heterogeneous slopes (`c.var`), `##`, or `newvar=` terms. Clustering is one-way only. Not every `reghdfe` solver or degrees-of-freedom option is available. |
| `civreghdfe` | `ivreghdfe` | aw/fw/pw | 2SLS, LIML, Fuller, k-class, two-step GMM, and CUE, with robust, one- or two-way clustered, and HAC standard errors. `absorb()` doesn't take heterogeneous slopes, and `center` isn't supported. A few estimator and VCE combinations aren't available; see the help file. |
| `cqreg` | `qreg` | — | IID, robust, or clustered standard errors. No weights and no `absorb()`; add fixed effects as factor variables instead (`i.firm`). `predict` supports `xb`, `residuals`, and `stdp`. It uses a different solver and density estimator than `qreg`, so coefficients and standard errors can differ slightly. It exits with an error if the solver doesn't converge. |
| `cpsmatch` | `psmatch2` | — | Estimates the ATT for one outcome. No ATE/ATU and no observation weights. Matching without replacement requires `neighbor(1)`. Standard errors use `psmatch2`'s default formula (below); the `ai()` option isn't supported. |
| `cppmlhdfe` | `ppmlhdfe` | aw/fw/pw | Supports fixed effects with interactions and heterogeneous slopes, exposure and offset, robust or multiway clustered standard errors, and `ppmlhdfe`'s separation checks (`fe`, `simplex`, `relu`, `mu`). Individual fixed effects can be saved, and `d()` saves their sum for `predict`. `ppmlhdfe`'s undocumented debugging options aren't available, and the iteration path can differ even when the final estimates agree. See [cppmlhdfe details](README_cppmlhdfe.md). |

## cpsmatch standard errors

`cpsmatch` reports the same ATT standard error as `psmatch2` does by default:

```
SE = sqrt(s1²/n1 + s0² · Σ w0² / n1²)
```

Here `n1` is the number of matched treated observations, `s1²` and `s0²` are the unweighted outcome variances among matched treated observations and the controls they're matched to, and `w0` are the controls' matching weights. Squaring the weights accounts for controls used more than once and for fractional matches. Like `psmatch2`'s default, this ignores the uncertainty from estimating the propensity score.

The Unmatched row also follows `psmatch2`: its standard error comes from an OLS regression of the outcome on treatment.

## String variables

Most commands handle numeric variables, str# variables up to `str2045`, and strL values of text up to 2,045 bytes. Longer or binary strL values aren't supported, and neither is creating strL variables. When a command can't handle a strL, it either passes the work to the Stata command it replaces (`cencode`, `cdecode`, `csplit`) or stops with an error (`csort`, `crangejoin`). `cimport` and `cexport` are the exception: they handle long text and binary data for the Excel, SAS, SPSS, and shapefile formats.

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

Speedups depend on the number of observations and variables, variable types, options, hardware, and thread count. In some cases ctools is slower than the command it replaces. See the [benchmarks](README_BENCHMARK.md) for measured results.

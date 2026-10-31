# chdfe

Remove one or more sets of fixed effects from numeric variables. `chdfe`
follows the user-written `hdfe` command, using the shared C solver that also
supports `creghdfe`. It returns transformed variables and absorbed degrees of
freedom; it does not fit a regression.

## Syntax

```stata
chdfe varlist [if] [in] [aw fw pw], absorb(absvars) generate(prefix) [sample(newvar) options]
chdfe varlist [if] [in] [aw fw pw], absorb(absvars) clear [keepvars(varlist) keepids options]
```

Use standard Stata weight syntax, such as `[aw=weight]`. Factor variables and
time-series operators are allowed in the variable list. Base and omitted
indicators are excluded from output.

`generate()` creates double variables by adding a prefix to each input name.
Expanded terms use the `hdfe` naming convention: names start with `__`, dots
become `__`, and `#` becomes `_X_`. The original data and row order are
preserved. `sample()` saves the zero-one transformation sample.

`clear` keeps the selected rows and transformed variables. It also retains the
source variables needed for expanded terms and simple weight variables. `keepvars()`
retains additional variables; `keepids` retains category IDs as `__ID1__`,
`__ID2__`, etc., and cluster IDs as `__CL1__`, etc. Failed calls restore the
original data.

## Fixed effects and sample

`absorb()` accepts categorical variables, `i.var`, and interactions such as
`firm#year`. It also accepts group-specific slopes: `firm#c.time` removes
slopes, while `firm##c.time` removes intercepts and slopes.
`firm##c.(time age)` allows several slopes. Categorical keys may be strings;
`strL` keys are limited to 2,045 bytes per value.
Fixed-effect coefficient output is not available.

Rows with missing inputs or categories, or missing or zero weights, are
excluded. Negative weights are rejected. Singleton rows are removed
repeatedly unless `keepsingletons` is specified. This follows `hdfe`, including
its removal of singleton rows with frequency weights greater than one.
Residuals are missing outside the retained sample. No retained rows produces
error 2000.

## Options

| Option | Meaning |
|---|---|
| `clustervars(varlist)` | Cluster categories or interactions for absorbed degrees of freedom |
| `dofadjustments(list)` | `all` (default), `none`, or a selection of `pairwise`/`firstpair`, `clusters`, and `continuous` |
| `groupvar(newvar)` | Connected-component IDs for the first pair of nonnested categorical effects |
| `tolerance(#)` | Convergence tolerance; default `1e-8` |
| `maxiterations(#)` | Iteration limit; default 10,000; `.` removes the practical limit |
| `poolsize(#)` | Columns projected together; default 10 |
| `transform(str)` | `symmetric_kaczmarz` (default), `kaczmarz`, or `cimmino` |
| `acceleration(str)` | `conjugate_gradient`/`cg` (default), `steepest_descent`/`sd`, `aitken`, or `none` |
| `threads(#)` | Maximum workers; zero uses the ctools setting |
| `verbose`, `verbose(#)`, `timeit` | Report time by computation stage; `verbose(0)` is quiet |

Conjugate gradient requires a symmetric transform. Iteration paths and counts
may differ from `hdfe`. Nonconvergence produces error 430. `groupvar()`
requires pairwise or first-pair adjustments; no output is created without an
eligible pair. With more than two categorical effects, the degrees-of-freedom
calculation is generally an upper bound, as in `hdfe`.

## Saved results

The reference results are `e(df_a)`, `e(N_hdfe)`, `e(N_hdfe_extended)`,
`e(mobility)`, `e(M_due_to_nested)`, and the per-effect scalars `e(Kj)`,
`e(Mj)`, `e(Mj_exact)`, `e(Mj_nested)`, and `e(Gj)`. `e(absvars)` and
`e(extended_absvars)` name the effects. `e(cmd)` is `chdfe`.

Additional results are `e(N)`, `e(num_singletons)`, `e(iteration)`,
`e(converged)`, and `e(dofmethod)`. `chdfe, version` reports the ctools release
and stores `e(version)`.

## Examples

```stata
chdfe sales capital, absorb(firm year) generate(within_)
chdfe sales capital [aw=weight], absorb(firm year) generate(r_) sample(used)
chdfe sales capital, absorb(firm##c.year) clear keepvars(firm year)
```

See also [creghdfe](creghdfe.md) and `help chdfe` in Stata.

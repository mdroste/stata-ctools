# creghdfe

Linear regression with high-dimensional fixed effects.

## Overview

`creghdfe` estimates linear regressions with multiple high-dimensional fixed effects. It removes fixed effects in C using a conjugate-gradient solver. See [Compatibility](../COMPATIBILITY.md) for differences from `reghdfe`.

## Syntax

```stata
creghdfe depvar indepvars [if] [in] [weight] [, absorb(varlist) options]
```

## Options

### Model Options
| Option | Description |
|--------|-------------|
| `absorb(absvars)` | Categorical fixed effects; omission estimates a model with an intercept and no absorbed effects |
| `absorb(id##c.x)`, `absorb(id#c.x)`, `absorb(id##c.(x z))`, `absorb(id##(c.x c.z))` | Heterogeneous slopes, with or without the intercept of `id`; with slopes alone the model has no constant, as in reghdfe |
| `absorb(newvar=absvar)` | Save that effect in `newvar` (its slopes in `newvarSlope1`, ...) |
| `absorb(absvars, savefe)` | Save the other fixed effects in `__hdfe#__` (slopes in `__hdfe#__Slope1`, ...); `generate` is a synonym. Every run first drops `__hdfe*`, as reghdfe, so reruns and `bootstrap:` work |
| `absorb(newvar=absvar, replace)`, `absorb(absvars, noisily)` | Replace an existing `newvar` target (`newvarSlope#` must be new); show the parsed effects. Other suboptions are r(198) |
| `absorb(ib2.var)`, `absorb(ibn.var)` | As `i.var`; `e(absvars)` keeps the term as typed |
| `residuals(newvar)` (`res()`, `resid()`, `resid2()`) | Save residuals in a new variable |
| `residuals` (`res`, `resid`) | Save residuals in `_reghdfe_resid`, replacing it, as reghdfe |
| `groupvar(newvar)` | Save mobility-group identifiers |
| `group(groupvar)` | One regression row per group (a patent); its variables must be constant within the group |
| `individual(indvar)`, `aggregation(mean\|sum)` | With `group()`: individual (inventor) fixed effects, `indvar` in `absorb()`, entering each group as the mean (default) or sum of its individuals' effects, as reghdfe |
| `dofadjustments(string)` | Absorbed degrees-of-freedom adjustments; default `all` |
| `tolerance(#)` | Projection tolerance; default `1e-8` |
| `iterate(#)` | Maximum projection iterations; default 16,000 |
| `transform(str)`, `acceleration(str)` | reghdfe's MAP variants with its names and abbreviations: `symmetric_kaczmarz` (default), `kaczmarz`, `cimmino`; `cg` (default), `sd`, `aitken`, `none`. `kaczmarz` with `cg` runs with reghdfe's warning; `rand_kaczmarz` is r(3498), as there |
| `technique(map\|lsmr\|lsqr)`, `preconditioner()` | Checked as reghdfe checks them (`technique(cg)` r(3498)); creghdfe always solves with MAP to the same tolerance |
| `prune`, `poolsize(#)`, `compact`, `fastregress`, `parallel(#)` | Accepted as reghdfe accepts them; no effect (data are held in C; threads replace worker processes) |
| `abort(0)` | reghdfe's undocumented solver property: keep the projection after `iterate()` iterations without convergence, `e(converged)` = 0. `miniter()`, `min_ok()`, `accel_start()`, `accel_freq()` are accepted with no effect; other `key(value)` options are r(198) |
| `version(#)` | Not supported: `version(3)` and `version(5)` (historical reghdfe) are r(198), other numbers r(9) as in reghdfe |
| `nopartialout`, `noregress`, `keepmata`, `varlist_is_touse` | Not supported (r(198)): they hand reghdfe's Mata `HDFE` object to other commands, and creghdfe has none |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |
| `snapshot`, `fast` | Use a recovery snapshot, or skip recovery; mutually exclusive |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce(vcetype)` | Variance-covariance estimator |

Where `vcetype` can be:
- `unadjusted` or `ols` - Standard VCE assuming homoskedasticity
- `robust` - Heteroskedasticity-robust (HC1) standard errors
- `cluster clustvars` - Cluster-robust standard errors; with several variables (up to ten), multi-way clustering as in reghdfe (Cameron-Gelbach-Miller inclusion-exclusion over every intersection, one small-sample factor with the fewest clusters, negative eigenvalues set to zero). When that correction leaves e(V) indefinite and `ereturn post` would set it to zero (reghdfe gives zero then, depending on the order of the cluster variables), creghdfe posts a positive semidefinite e(V) with the same coefficient block; only the constant's row and column differ
- `dkraay [#]` - Driscoll-Kraay standard errors over the `tsset`/`xtset` time variable, bandwidth `#` (lags + 1; reghdfe's default when omitted), `e(dkraay_bw)`
- `bootstrap [, bootopts]` - Bootstrap standard errors, computed in the plugin (reghdfe has no such option)

`vce(bootstrap)` resamples observations, or clusters with `cluster()`, within
`strata()` if given, and posts what `vce(bootstrap)` posts after official
commands, so `estat bootstrap` works. Inference is z-based (no `e(df_r)`); the
fit's own covariance is `e(V_modelbased)`. Its options are `reps(#)` (default
50), `seed(#)`, `cluster(clustvar)`, `strata(varlist)`, `draws(stata|fast)` (default `stata`),
`mse`, `level(#)`, `nodots`/`dots(#)` and
`saving(filename [, replace double every(#)])`. As with `bootstrap`, `mse`
posts `e(vcetype)` `Bstrap *` and `e(mse)`, and the command's own `level()`
sets the bootstrap's (the two may not differ); the table at estimation shows
that level.

`saving()` writes the replicates' coefficients as `bootstrap` writes them: one
observation per replicate, one variable per term of `e(b)` (`_b_name`, or
`_bs_#`), with `bootstrap`'s labels and characteristics, so `bstat` reads the
file. Failed replicates are missing throughout, omitted terms zero. Values
are floats unless `double` is given; `e(V)` uses the double replicates either
way.

Random numbers come from Stata's generator: `set seed`, `seed()`, `set rng` and
`set rngstream` apply, and `e(rngstate)` records the state. With
`draws(stata)`, the default, the draws are exactly `bsample`'s, so the results
reproduce `bootstrap _b, reps(#) seed(#): reghdfe ...` for the same seed.
`draws(fast)` takes a master seed from Stata's generator and draws in the
plugin: reproducible from the seed and faster for large observation-level
bootstraps, but not `bsample`'s draws. Weights are allowed; without `cluster()`, fweights are resampled as
the observations they stand for. With `draws(stata)` those are `bsample`'s
draws on the fweight-expanded data (in `expand`'s order), so the results are
the prefix's after `expand`. Within `strata()`, with observations resampled,
`draws(stata)` also takes the uniforms the prefix's `bsample` takes after its
draws (one per observation, without changing the sample), so that later
replicates and the final RNG state match the prefix. Heterogeneous slopes are
re-estimated in each replicate, as the fit estimates them.

### Fixed-Effect Standard Errors
| Option | Description |
|--------|-------------|
| `fevce(bootstrap [, feopts])` | Bootstrap standard errors of the `savefe` fixed effects |

`fevce(bootstrap)` requires `absorb(..., savefe)`. Each replicate re-estimates
the model and recovers the effects of the levels in its resample; a level's
standard error is the standard deviation of its effect over the replicates that
estimated it, saved in `__fese1__`, `__fese2__`, ... on the estimation sample
(missing where fewer than two replicates did). With `vce(bootstrap)` both share
one set of replicates, resampling options given to both must agree, and
`level()`, `mse` and `dots()` are those of `vce()` (`fevce()` may repeat them). Alone,
`fevce()` leaves `e(V)` as `vce()` computes it. `feopts` are the `vce(bootstrap)`
options (`saving()` goes in `vce(bootstrap)` when both are given) and:

- `sestub(stub)` - name the standard errors `stub1`, `stub2`, ...
- `nreps` - also save each level's count of estimating replicates (`__fenrep#__`, or `stub#_nrep`)
- `ci` - also save percentile bounds at `level()` (`__felb#__` and `__feub#__`, or `stub#_lb` and `stub#_ub`; `e(fe_lb)`, `e(fe_ub)`, `e(fe_level)`). Every replicate's effects are kept (4 bytes each); r(909) when they do not fit in memory
- `dse(newvar)` - save the standard error of the summed fixed effects d

With heterogeneous slopes, a slope saved in `__hdfe#__Slope1` has its standard
error in `__fese#__Slope1` (and `__fenrep#__Slope1`, `__felb#__Slope1`, ...). A
slope is estimable in a level where the replicate's slope variable varies; the
connected sets and the normalization concern the intercepts.

A level is estimable in a replicate when it was drawn and its full-sample
connected set did not split; replicates whose omitted regressors differ from
the full sample's estimate no fixed effects. `e(fe_bs_table)` summarizes, per
effect, the levels never and always estimable and the minimum and median
replicate counts.

With two or more absorbed effects, individual fixed effects are identified only
up to a normalization within each connected set, so their bootstrap standard
errors depend on it. `creghdfe` displays a warning (also on replay) and sets
`e(fe_norm_warning)` to 1. The summed effect d is identified: `dse()` gives
standard errors that do not depend on the normalization.

### Wild Bootstrap
| Option | Description |
|--------|-------------|
| `wildbootstrap[(wildopts)]` | Wild (cluster) bootstrap p-values and confidence intervals of the coefficients |

`wildbootstrap` keeps e(b), e(V) and the analytic table and adds a table of
bootstrap p-values (H0: b = 0) and test-inversion intervals. It needs
`vce(robust)` (one weight per observation) or `vce(cluster clustvar)`. The
null is imposed by default (WCR). The conventions are boottest's: for the
same `seed()` the weights and the p-values are boottest's. Options:
`reps(#)` (999), `seed(#)`, `weighttype(rademacher|mammen|webb|normal|gamma)`,
`nonull` (WCU), `coefs(names)`, `ptype(symmetric|equaltail|lower|upper)`,
`noenumerate`, `fedgp(full|partialled)`, `noci`, `level(#)`,
`bootcluster(varname)` and `nodots`/`dots(#)`.

- With Rademacher weights and 2^C <= reps, all 2^C sign vectors are used
  (e(wb_enumerated) = 1; no random draws).
- Effects nested in the bootstrap clusters cost a few cluster sums per
  replication. An effect that crosses them is partialled out of each
  replicate (the full-model bootstrap, as boottest after reghdfe);
  `fedgp(partialled)` follows Stata's `wildbootstrap areg` instead.
- Intervals are bounded exactly where the p-value crosses the level;
  boottest stops within one replication's step of it.
- Not with `vce(bootstrap)`; with `vce(robust)` and fweights it needs
  `bootcluster()`.

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information and timing |
| `nowildtable` | Suppress the wild bootstrap table |
| `noheader` (`nohe`), `notable`, `nofootnote` | Suppress the header, the table, the absorbed degrees of freedom |
| `nosample` | Do not save `e(sample)` |
| `timeit` | reghdfe's timer tables with creghdfe's steps; `r(timer)` |
| `nowarn` | Omit the `keepsingletons` warning |

`vce()` reads its type as reghdfe does: a word starting with `r` is `robust`, with `un` or `ols`
`unadjusted` (as is `conventional`), with `cl` `cluster` and with `dk` `dkraay`; any other is r(9).
The header and the absorbed-degrees-of-freedom table have reghdfe's layout, and its notes and
warnings come in its order. A factor-variable dependent variable naming one indicator
(`creghdfe 1.foreign x, absorb(id)`) is allowed, as in reghdfe; `i.foreign` is r(9).

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

* Two-way clustering; Driscoll-Kraay standard errors
creghdfe y x1 x2, absorb(firm year) vce(cluster firm year)
xtset firm year
creghdfe y x1 x2, absorb(firm) vce(dkraay 3)

* Firm intercepts and trends, saved as fe and feSlope1
creghdfe y x1 x2, absorb(fe=firm##c.trend year)

* With robust standard errors
creghdfe y x1 x2, absorb(firm year) vce(robust)

* With cluster bootstrap standard errors, reproducing the bootstrap prefix
creghdfe y x1 x2, absorb(firm year) vce(bootstrap, reps(500) seed(1) cluster(firm) draws(stata))
estat bootstrap, percentile

* Wild cluster bootstrap p-values and intervals (boottest's conventions)
creghdfe y x1 x2, absorb(firm year) vce(cluster state) wildbootstrap(reps(9999) seed(1) weighttype(webb))

* Bootstrap standard errors of saved firm effects, sharing the replicates
creghdfe y x1 x2, absorb(firm, savefe) vce(bootstrap, reps(200) seed(1)) fevce(bootstrap, nreps)

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

### Implementation details

**CG Solver and HDFE Absorption:**
- **Symmetric Kaczmarz transform**: Forward and backward projection sweeps provide the symmetric operator used by CG
- **Pre-computed inverse group counts**: Multiplies by `1.0/count` instead of dividing, to avoid repeated division
- **Counting sort for CSR index construction**: Builds compressed sparse row indices in O(N+L) instead of O(N log N), where L is the number of FE levels
- **Queue-based singleton removal**: Counts rows per FE level once, then peels singleton levels from a queue until none remain, so chains where dropping one singleton creates another cost O(N·G) in total

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
4. **VCE Computation**: Calculate standard errors using the selected covariance estimator and degrees-of-freedom adjustments

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
| `e(N_full)` | Observations before singletons were dropped |
| `e(converged)`, `e(ic)` | 1 when the projection converged (0 after `abort(0)` without convergence, with `e(ic)` 0, as reghdfe's MAP); iterations |
| `e(N_clust)` | Fewest clusters of any cluster variable (if clustered; time periods with `vce(dkraay)`) |
| `e(N_clustervars)`, `e(N_clust#)` | Number of cluster variables, and the clusters of each |
| `e(dkraay_bw)` | Driscoll-Kraay bandwidth (if `vce(dkraay)`) |
| `e(N_hdfe)`, `e(N_hdfe_extended)` | Absorbed effects, and their intercepts and slopes |
| `e(fe_N_reps)`, `e(fe_split_reps)` | Replicates run, and replicates in which a connected set split (if `fevce(bootstrap)`) |
| `e(fe_norm_warning)` | 1 with two or more absorbed effects (if `fevce(bootstrap)`) |
| `e(wb_reps)`, `e(wb_misreps)`, `e(wb_enumerated)`, `e(wb_N_clust)`, `e(wb_level)` | Wild bootstrap replications, infeasible ones, enumeration, clusters, level |

### Macros
| Result | Description |
|--------|-------------|
| `e(cmd)` | `creghdfe` |
| `e(depvar)` | Dependent variable name |
| `e(indepvars)` | Regressors partialled out, as reghdfe lists them (`2bn.rep78`, ...; base levels left out), and `_cons` when reported |
| `e(estat_cmd)`, `e(footnote)` | `_creghdfe_estat` (`estat summarize`, `vce`, `ic` with `e(df_m)` + 1 degrees of freedom, `bootstrap`) and `_creghdfe_footnote` (the absorbed degrees of freedom, for prefix displays) |
| `e(absorb)` | Absorbed fixed effects |
| `e(vce)` | VCE type |
| `e(clustvar)`, `e(clustvar#)` | Cluster variables (if clustered; the time variable with `vce(dkraay)`) |
| `e(absvars)`, `e(extended_absvars)` | Absorbed effects as specified, and their intercepts and slopes (`id`, `id#c.x`) |
| `e(title3)`, `e(title4)` | Robustness titles of the header (clustered and Driscoll-Kraay) |
| `e(fevce)` | `bootstrap` (if `fevce(bootstrap)`) |
| `e(fe_se)`, `e(fe_nrep)`, `e(fe_dse)` | Variables with the FE standard errors, replicate counts, and d's standard errors |
| `e(fe_rngstate)`, `e(fe_draws)` | RNG state before the first draw, and `draws()` |
| `e(wb_method)`, `e(wb_weighttype)`, `e(wb_ptype)`, `e(wb_coefs)`, `e(wb_bootclust)`, `e(wb_fe_path)`, `e(wb_rngstate)` | Wild bootstrap settings |

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix |
| `e(dof_table)` | Per absorbed intercept and slope: categories, redundant, coefficients, inexact, nested; rows `1.id`, `1.a#1.b`, `1.id#c.x` and columns `Categories`, `Redundant`, `Num Coefs`, `Inexact?`, `Nested?`, as reghdfe names them |
| `e(wb_t)`, `e(wb_p)`, `e(wb_ci)`, `e(wb_cistatus)` | Analytic t, wild bootstrap p-values, intervals and their status (open lower, open upper, disjoint, empty, failed) |
| `e(fe_bs_table)` | Per effect: levels, never/always estimable, min and median replicate counts (if `fevce(bootstrap)`) |

## Technical Notes

- The absorbed fixed effects consume degrees of freedom
- Singleton groups (groups with only one observation) are dropped, repeatedly, until none remain; with frequency weights a group is a singleton only when its total weight is 1, as in reghdfe. Factor variables are expanded on the sample after singleton removal, as in reghdfe, so levels and base levels follow the estimation sample
- The absorbed degrees of freedom account for redundant fixed-effect coefficients according to `dofadjustments()`
- The constant term is absorbed with the fixed effects

## Comparison with reghdfe

| Feature | `reghdfe` | `creghdfe` |
|---------|-----------|------------|
| Implementation | Stata + Mata | C with OpenMP |
| Algorithm | Method of Alternating Projections | CG with Kaczmarz |
| Parallelization | Limited | Yes |

Memory use depends on the dataset, options, and number of threads.

## Degrees of Freedom

The residual degrees of freedom is computed as:
```
df_r = N - K - e(df_a)
```
for an unclustered model, where K is the number of retained slope coefficients and `e(df_a)` already accounts for redundant fixed effects. Clustered inference uses the cluster degrees of freedom.

Redundant coefficients are counted as `reghdfe` counts them. Effects nested within any cluster variable (the time variable with `vce(dkraay)`) cost no degrees of freedom and are set aside. Each remaining effect after the first then loses the largest number of mobility groups it forms with any earlier remaining effect, and at least one coefficient. `dofadjustments()` takes `reghdfe`'s values: `all` (default), `none`, or any of `firstpair` or `pairwise`, `clusters` and `continuous`. With `continuous`, a heterogeneous slope is redundant in each level where its standardized variable is constant (or zero, without an intercept). Slopes are partialled out as separate factors (each centered within its levels when its effect has an intercept), which converge to the same projection as reghdfe's joint fit; their saved estimates are recovered with reghdfe's joint fit per level. Slopes work with `vce(bootstrap)`, `fevce()` and `wildbootstrap()`. Group and individual fixed effects (`group()`, `individual()`, `aggregation()`) follow reghdfe: one row per group, one coefficient per individual with none counted as redundant, and groups with an individual found in no other group dropped as singletons. The projection is solved with LSMR (reghdfe's `technique(lsmr)`), with diagonal preconditioning, 100 times more tightly than `tolerance()`. With `individual()`, fweights, slopes, saved effects, `vce(bootstrap)`, `fevce()`, `wildbootstrap()` and `groupvar()` are refused. With aweights or pweights the estimates are the weighted least-squares fit, equal to the regression on explicit individual-share variables; reghdfe 6.13.1 does not reach it there. creghdfe exits with r(459) when a group and an individual do not identify the observations, as reghdfe means to (its check never runs).

## See Also

- [ctools Overview](../../README.md)
- [reghdfe](http://scorreia.com/software/reghdfe/) - Original Stata implementation
- [cqreg](cqreg.md) - Quantile regression
- [civreghdfe](civreghdfe.md) - IV regression with HDFE
- [HDFE benchmark datasets](../performance/HDFE_DATASETS.md) - `creghdfe` compared with `reghdfe` on 15 benchmark graphs

## Validation and failure behavior

residuals(), groupvar(), and savefe output names must be distinct new variables, with reghdfe's exceptions: every run drops `__hdfe*`, `residuals` without a name replaces `_reghdfe_resid`, and `absorb(..., replace)` replaces existing `newvar=` targets (after the other names are checked, so a failed check replaces nothing). Other existing variables, including inputs and fixed-effect identifiers, are never replaced. Failed or interrupted operations restore the dataset without writing a snapshot to disk (the variables creghdfe added are dropped again); option `snapshot` uses `preserve` instead and option `fast` undoes nothing. iterate() and tolerance() must be positive. If the fixed-effect solver reaches its iteration limit before converging, the command returns error 430 and does not store new estimates.

# civreghdfe

Instrumental variables regression with high-dimensional fixed effects.

## Overview

`civreghdfe` estimates instrumental variables models in C, with optional absorbed fixed effects. It defaults to two-stage least squares (2SLS) and also supports LIML, Fuller LIML, k-class, two-step GMM, and continuously updated GMM.

## Syntax

```stata
civreghdfe depvar (endogvars = instruments) [exogvars] [if] [in] [weight] [, absorb(varlist) options]
civreghdfe [, first ffirst rf sfirst dropfirst droprf dropsfirst level(#) noheader nofooter eform(string) plus]
```

Without a variable list (or with `estimates replay`), `civreghdfe` redisplays the last results, with `ivreg2`'s replay options: `first`, `rf` and `sfirst` show the equations stored by `savefirst`, `saverf` and `savesfirst`, `ffirst` the summary of `e(first)`, and `drop*` drops the stored equations. `by:` with replay is error 601; replay after another command is error 301.

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `absorb(varlist)` | Optional categorical fixed effects; heterogeneous slopes as reghdfe absorbs them (`id#c.x`, `id##c.x`, `id##c.(x z)`, at most ten intercepts and slopes), not with `vce(bootstrap)`, `fevce()` or `wildbootstrap()` |
| `absorb(newvar=absvar ...)`, `absorb(..., savefe)` | Save the estimated effects (in `newvar`, or `__hdfe#__`), as reghdfe; the effects of y - Xb - e, normalized as reghdfe's with the constant in the first effect as ivreghdfe saves them; with `partial()`, the partialled regressors' coefficients are recovered by FWL, so the effects are those of the model without it |

### Variance/Standard Error Options
| Option | Description |
|--------|-------------|
| `vce(vcetype)` | Variance estimation: `unadjusted`, `robust`, `cluster clustvar`, or `cluster clustvar1 clustvar2` |
| `robust`, `cluster(varlist)` | Aliases for `vce(robust)` and `vce(cluster varlist)` |
| `bw(#)` | Kernel bandwidth: AC standard errors, or HAC with `robust` (data must be `tsset`); `bw(auto)` selects it (Newey-West 1994, time series only) |
| `kernel(string)` | Kernel for `bw()`: `bartlett` (default), `parzen`, `qs`, `truncated`, `thann`, `thamm` (Tukey-Hamming), `dan` (Daniell), `ten` (Tent) |
| `center` | Centered moments in the robust, cluster and kernel S matrix |
| `dkraay(#)` | Driscoll-Kraay standard errors with bandwidth `#` (tsset panel data) |
| `kiefer` | Kiefer standard errors: within-panel autocorrelation (tsset panel data) |
| `vce(bootstrap [, bootopts])` | Pairs bootstrap standard errors, run in the plugin (below) |
| `fevce(bootstrap [, feopts])` | Bootstrap standard errors of the saved fixed effects (below) |
| `wildbootstrap[(wildopts)]` | Wild restricted efficient (WRE) bootstrap p-values and confidence intervals (below) |
| `nowildtable` | Suppress the wild bootstrap table |

`vce(bootstrap)` resamples observations, or clusters with `cluster()`, within
`strata()` if given, and posts what `vce(bootstrap)` posts after official
commands (`estat bootstrap` works; z-based inference, no `e(df_r)`, a Wald
chi2 test in the header). Each replicate re-estimates the model with the
chosen estimator after dropping its own singletons, partialling out the
fixed effects and `partial()` again and applying the fit's collinearity rules;
a replicate whose omitted regressors or instruments differ from the full
sample's is missing. The `bootopts` are `creghdfe`'s: `reps()`, `seed()`,
`cluster()`, `strata()`, `draws(stata|fast)` (default `stata`), `mse`, `level()`, `nodots`,
`dots()`; `mse` and `level()` behave as in creghdfe.

The fit keeps a base VCE for what needs one: the identification and
overidentification statistics, the first-stage and reduced-form regressions,
and the moment covariance S of `gmm2s` and `cue` use the cluster-robust VCE on
the bootstrap's `cluster()`, otherwise the robust one; it is posted as
`e(V_modelbased)`. In `gmm2s` and `cue` replicates, S is that of the data the
replicate stands for: each drawn copy of a cluster is its own cluster (as
`idcluster()` makes it). With `draws(stata)` and the same seed the results equal
the bootstrap of `ivreghdfe` on the drawn, expanded data (the validation's
oracle; `ivreghdfe` cannot run under the `bootstrap` prefix itself). Weights are
allowed: aweights and pweights multiply the resampled rows; without
`cluster()`, fweights are resampled as the observations they stand for (with
`draws(stata)`, `bsample`'s draws on the expanded data). `saving()` saves the
replicates as `bootstrap` does. Not with `bw()`, `kernel()`, `dkraay()`, `kiefer` (resample panels
with `cluster()`), `b0()`, `smatrix()` or `wmatrix()`.

`fevce(bootstrap)` needs saved effects (`savefe` or `newvar=`): each replicate
recovers the effects of y - Xb - e as the fit saves them, and the standard
error of `__hdfe#__` is saved in `__fese#__` (of `newvar` in `newvar_se`);
`sestub()`, `nreps`, `ci` and `dse()` as in `creghdfe`, also with
`partial()`. With `vce(bootstrap)` both
share one set of replicates (`level()`, `mse` and `dots()` are those of
`vce()`); alone, `e(V)` is that of `vce()`, which with
`gmm2s` or `cue` must be robust or one-way cluster. With two or more absorbed
effects a warning notes that individual effects depend on the normalization;
their sum (`dse()`) does not.

`wildbootstrap` adds bootstrap p-values of H0: b = 0 and confidence intervals
by test inversion for the coefficients (endogenous and exogenous), after the
analytic fit, which it leaves unchanged. It is the wild restricted efficient
bootstrap of Davidson and MacKinnon (2010) as `boottest` computes it: the
restricted structural equation is estimated by LIML, the first stage
regresses the endogenous regressors on the instruments and that residual,
and each replication perturbs both residuals by the same weight, then
re-estimates with the command's estimator (2SLS, LIML, Fuller, k-class) and
its robust or cluster-robust standard error. `nonull` bootstraps from the
unrestricted LIML fit (WUE). The `wildopts` are `creghdfe`'s; with the same
`seed()` the p-values equal `boottest`'s after `ivreg2` with the absorbed
effects as dummies (`boottest` after `ivreghdfe` with `absorb()` does not
reproduce the WRE statistic, and it does not take `fuller()` or `kclass()`).
Effects nested in the bootstrap clusters cost nothing; others are partialled
out of each replicate's perturbed data. `partial()` variables stay regressors
of the bootstrap model. Needs `vce(robust)` or `vce(cluster)`; not with
`gmm2s`, `cue`, kernel VCEs, `b0()`, `smatrix()`, `wmatrix()` or
`vce(bootstrap)`; with `fevce(bootstrap)` both run on one fit.

### Estimation Options
| Option | Description |
|--------|-------------|
| `tolerance(#)` | CG solver convergence tolerance (default: 1e-8) |
| `maxiter(#)` | Maximum CG iterations (default: 500) |
| `liml`, `fuller(#)`, `kclass(#)` | Alternative IV estimators |
| `gmm2s`, `cue` | Two-step or continuously updated GMM |
| `wmatrix(matrix)`, `smatrix(matrix)` | `ivreg2`'s user-supplied weighting matrix W (the `gmmw` model; first step of `gmm2s`) and covariance matrix S |
| `partial(varlist)` | Partial out specified model terms; `_cons` partials out the constant (also alone; not with `absorb()` or `noconstant`), `_all` every included exogenous regressor; implies `bvclean`, as in `ivreg2`; `fwl()` is an alias |
| `psd0`, `psda` | Make the moment covariance S positive semidefinite (negative eigenvalues set to 0, or to their absolute values), as `ivreg2`; not with the bootstraps |
| `tvar()`, `ivar()` | `ivreg2`'s time and panel variables: they must be the `tsset` ones (error 5 otherwise, as `ivreg2`) |
| `bvclean` | `e(b)` and `e(V)` without base levels and omitted terms (factor levels named as `ivreg2` names them, e.g. `1bn.a`); not with the bootstraps |
| `fvall`, `fvsep`, `useqr` | Accepted for `ivreg2` compatibility; they do not change civreghdfe's results |
| `groupvar(newvar)` | `reghdfe`'s mobility group of the first two absorbed effects on `e(sample)` (with two or more effects) |
| `technique()`, `transform()`, `acceleration()`, `preconditioner()`, `prune`, `poolsize()`, `nosample`, `nowarn` | `reghdfe`'s options that `ivreghdfe` passes on: the names are checked as `reghdfe` checks them; civreghdfe's solver does not use them, and `nosample` keeps `e(sample)` (as in `ivreghdfe`) |
| `orthog()`, `endogtest()`, `redundant()` | Instrument orthogonality (excluded or included instruments), regressor endogeneity, and instrument redundancy tests |
| `residuals(newvar)` | Save residuals |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |

### Reporting Options
| Option | Description |
|--------|-------------|
| `first` | Report the first-stage regressions and `ivreg2`'s first-stage statistics |
| `ffirst` | Report `ivreg2`'s first-stage summary: F, Sanderson-Windmeijer, identification and Anderson-Rubin/Stock-Wright statistics |
| `sfirst` | Report the first-stage and reduced-form regressions as one system; `savesfirst`, `savesfprefix()` store it |
| `savefirst`, `saverf` | Store the first-stage and reduced-form regressions |
| `verbose`, `verbose(#)` | Display progress information (`verbose(#)` with # > 0) |
| `timeit` | Display the elapsed times |
| `small` (`sm`) | Use small-sample adjustments; `nohe`/`nofo` abbreviate `noheader`/`nofooter`, as in `ivreg2` |

### Weights
Supports `aweight`, `fweight`, and `pweight`, and `iweight` as `ivreg2` does: without `absorb()` (error 101, as `ivreghdfe`), with the homoskedastic VCE only (error 101 with `robust`, `cluster()`, `gmm2s` or kernels), and not negative (error 402). `e(N)` and the residual degrees of freedom use the integer part of their sum; with a non-integer sum, `ivreg2` also truncates N in its identification statistics, which civreghdfe computes with the exact sum.

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

* Kernel-based standard errors (the data must be tsset; nlswork is xtset idcode year)
civreghdfe ln_wage (tenure = union) age, absorb(idcode) bw(2)            // AC
civreghdfe ln_wage (tenure = union) age, absorb(idcode) bw(2) robust     // HAC
civreghdfe ln_wage (tenure = union) age, absorb(idcode) dkraay(2)        // Driscoll-Kraay
civreghdfe ln_wage (tenure = union) age, absorb(idcode) kiefer

* Display first-stage statistics
civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) first

* With weights
civreghdfe ln_wage (tenure = union) age [aw=hours], absorb(idcode)

* Multiple fixed effects
civreghdfe y (x1 = z1 z2) x2 x3, absorb(firm year industry)

* Verbose timing output
civreghdfe ln_wage (tenure = union) age, absorb(idcode) verbose
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
Compute residuals using the original endogenous regressors after any fixed-effect removal. Calculate the covariance matrix using the selected estimator and variance option.

## First-Stage Statistics

With the `first` option, `civreghdfe` reports:
- First-stage F-statistics for each endogenous variable
- Reports excluded-instrument relevance; interpretation depends on the model and covariance estimator

With the unadjusted covariance this is the classical F from the partial R-squared of the
excluded instruments. With robust, cluster, kernel (HAC/AC, Driscoll-Kraay, Kiefer) or
probability-weighted fits it is the Wald test of the excluded instruments under that
covariance, scaled as `ivreg2`'s `e(first)["F",]`; with clusters its denominator degrees of
freedom are G-1.

## Variance-Covariance Estimation

| VCE Type | Description |
|----------|-------------|
| `unadjusted` | Standard VCE assuming homoskedasticity |
| `robust` | Heteroskedasticity-robust (HC1) |
| `cluster clustvar` | Cluster-robust, clustered on `clustvar` |
| `cluster clustvar1 clustvar2` | Two-way cluster-robust (Cameron-Gelbach-Miller; small-sample factor from the smaller dimension) |
| `bw(#)` | AC: autocorrelation-consistent, conditionally homoskedastic |
| `bw(#) robust` | HAC: heteroskedasticity- and autocorrelation-consistent |
| `dkraay(#)` or `cluster(timevar) bw(#)` | Driscoll-Kraay: time-period clusters with a kernel over periods |
| `cluster(panelvar timevar) bw(#)` | Two-way clustering on panel and time with the kernel over time |
| `kiefer` | Truncated kernel over the whole panel (AC; HAC with pweights) |

Covariance calculations use the selected estimator and variance option. Residuals use the original endogenous regressors after any fixed-effect removal.

The kernel-based estimators follow `ivreg2`: the data must be `tsset` (error 111 otherwise; `dkraay` and `kiefer` need panel and time variables, error 5); kernel lags are differences of the time variable within panels, so gaps and unbalanced panels are handled like `ivreg2`'s lag operators and the panel structure does not depend on `absorb()`. `kernel()` requires `bw()` (error 102), fweights are not allowed with a kernel (error 101), a kernel with clustering must cluster on the time variable (or on panel and time), and the bandwidth may not exceed the time span of the sample. All of `ivreg2`'s kernels are available (the Daniell kernel is `e(kernel)` "Danielle", as in `ivreg2`; Quadratic Spectral, Daniell and Tent use every lag up to the time span). `bw(auto)` is Newey and West's (1994) selection from the first-step residuals, as `ivreg2`'s: Bartlett, Parzen and Quadratic Spectral kernels (bw=1 with a warning for others), time series only (error 198 with panel data), not with LIML, k-class or `b0()`; the choice is in `e(bwchoice)`. `center` centers the moments of the robust, cluster and kernel S as `ivreg2` does (the identification statistics are not centered, as in `ivreg2`).

The J statistic, the `orthog()` C statistic, the `endogtest()` statistic and the `gmm2s`/`cue` weighting matrix use the moment covariance S of the chosen VCE (robust, one- and two-way cluster, AC/HAC, Driscoll-Kraay, Kiefer). The Kleibergen-Paap rk LM and Wald statistics (`e(idstat)`, `e(widstat)`) are computed as `ranktest` does for any number of endogenous regressors. As in `ivreghdfe`, the identification and redundancy statistics and the model re-estimated for `endogtest()` use the Bartlett kernel whatever `kernel()` is, and are iid with `kiefer`. When S is rank deficient (e.g. fewer clusters than instruments), the J, C and endogeneity statistics are missing, as in `ivreg2`, and `gmm2s`/`cue` stop with error 506 (`ivreg2` stops as well: r(506) for `gmm2s`, an optimizer error for `cue`).

## Stored Results

`civreghdfe` stores the following in `e()`:

### Scalars
| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(df_r)` | Residual degrees of freedom (not with `vce(bootstrap)`) |
| `e(chi2)`, `e(p)` | Wald test of the regressors (with `vce(bootstrap)`, in place of `e(F)`) |
| `e(N_reps)`, `e(N_misreps)` | Bootstrap replicates, and those incomplete (with `vce(bootstrap)`) |
| `e(fe_N_reps)`, `e(fe_split_reps)`, `e(fe_norm_warning)` | As `creghdfe` (with `fevce(bootstrap)`) |
| `e(wb_reps)`, `e(wb_misreps)`, `e(wb_enumerated)`, `e(wb_N_clust)`, `e(wb_level)` | As `creghdfe` (with `wildbootstrap`) |
| `e(df_a)` | Absorbed degrees of freedom |
| `e(K)` | Number of regressors |
| `e(K_endog)` | Number of endogenous regressors |
| `e(K_exog)` | Number of exogenous regressors |
| `e(K_iv)` | Number of instruments (including exogenous) |
| `e(G)` | Number of absorbed intercepts and slopes |
| `e(N_hdfe)`, `e(N_hdfe_extended)` | Number of absorbed effects, and of their intercepts and slopes (as reghdfe) |
| `e(N_clust)` | Number of clusters (if clustered; time periods for Driscoll-Kraay) |
| `e(bw)` | Kernel bandwidth (kernel-based VCE; the time span for `kiefer`) |
| `e(F_first1)` | First-stage F-stat for 1st endogenous var |
| `e(F_first2)` | First-stage F-stat for 2nd endogenous var |
| `e(arf)`, `e(arfp)`, `e(ardf)`, `e(ardf_r)`, `e(archi2)`, `e(archi2p)` | Anderson-Rubin Wald tests (with `first`, `ffirst`, `sfirst`, `savefirst`, `savesfirst`) |
| `e(sstat)`, `e(sstatp)`, `e(sstatdf)` | Stock-Wright LM S statistic |
| `e(center)` | 1 with `center` |

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
| `e(extended_absvars)` | Absorbed intercepts and slopes, as reghdfe names them |
| `e(savefe)` | Variables holding the saved fixed effects |
| `e(vcetype)` | VCE type (`Bootstrap` with `vce(bootstrap)`) |
| `e(model)` | `ivreg2`'s estimator: `iv`, `ols` (no endogenous regressors), `liml` (also Fuller's), `kclass`, `gmm2s`, `cue` (also `b0()`), `gmmw` |
| `e(title)`, `e(subtitle)`, `e(hacsubtitleB)`, `e(hacsubtitleV)`, `e(hacsubtitleV2)` | `ivreg2`'s titles (`title()` sets `e(title)`) |
| `e(vce)` | `ivreg2`'s VCE description (`robust`, `cluster`, `two-way cluster`, `hac`/`ac` with kernel and bandwidth, `psd0`/`psda`) |
| `e(depvar0)`, `e(instd0)`, `e(inexog0)`, `e(exexog0)`, `e(insts0)` and the `1` lists, `e(partial1)` | `ivreg2`'s lists with base levels and omitted terms (`0`) and without them (`1`); partialled-out terms are not in them |
| `e(clist)`, `e(elist)`, `e(redlist)`, `e(partialsmall)` | `ivreg2`'s tested-variable lists and small-sample flag |
| `e(absvars)`, `e(title2)`, `e(dofmethod)` | `reghdfe`'s footnote results (with `absorb()`) |
| `e(vce)`, `e(prefix)`, `e(rngstate)`, `e(draws)`, `e(cluster_bs)`, `e(strata_bs)` | `bootstrap` and its settings (with `vce(bootstrap)`) |
| `e(fevce)`, `e(fe_se)`, `e(fe_nrep)`, `e(fe_dse)` | As `creghdfe` (with `fevce(bootstrap)`) |
| `e(wb_method)` | `wre`, or `wue` with `nonull` (with `wildbootstrap`) |
| `e(wb_weighttype)`, `e(wb_ptype)`, `e(wb_coefs)`, `e(wb_bootclust)`, `e(wb_fe_path)`, `e(wb_rngstate)` | As `creghdfe` (with `wildbootstrap`) |
| `e(clustvar)` | Cluster variable (if clustered) |
| `e(kernel)` | Kernel of a kernel-based VCE |
| `e(tvar)`, `e(ivar)` | tsset time and panel variables (kernel-based VCE) |
| `e(bwchoice)` | Bandwidth choice of `bw(auto)` |
| `e(firsteqs)`, `e(rfeq)`, `e(sfirsteq)` | Names of the stored first-stage, reduced-form and system equations (default prefix `_ivreg2_`, as `ivreg2`; each stored equation has `e(S)`, `ivreg2`'s subtitles, and `civreghdfe` replays it after `estimates restore`) |

### Matrices
| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix (the bootstrap's with `vce(bootstrap)`) |
| `e(V_modelbased)` | The base VCE's covariance (with `vce(bootstrap)`) |
| `e(b_bs)`, `e(reps)`, `e(bias)`, `e(se)`, `e(z0)`, `e(ci_normal)`, `e(ci_percentile)`, `e(ci_bc)`, `e(boot_status)` | As `creghdfe` (with `vce(bootstrap)`) |
| `e(fe_bs_table)` | As `creghdfe`, one row per saved effect (with `fevce(bootstrap)`) |
| `e(wb_t)`, `e(wb_p)`, `e(wb_ci)`, `e(wb_cistatus)` | As `creghdfe` (with `wildbootstrap`) |
| `e(first)` | `ivreg2`'s first-stage statistics per endogenous regressor (rmse, Shea and partial R2, F, Sanderson-Windmeijer and Angrist-Pischke statistics) |
| `e(S)`, `e(W)` | `ivreg2`'s covariance of the orthogonality conditions and weighting matrix (no `e(W)` for LIML and k-class) over the instruments kept, excluded first |
| `e(ccev)`, `e(cdev)` | Squared canonical correlations and Cragg-Donald eigenvalues |
| `e(dof_table)` | `reghdfe`'s absorbed degrees of freedom table (with `absorb()`) |

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
- Regressors and instruments absorbed by the fixed effects (sum of squares after absorption at most min(1e-6, `tolerance()`/10) times the sum of squared deviations before it) are omitted. `ivreghdfe` checks collinearity before absorbing, so it keeps such a column as partialling noise in its ranks unless its solver removes it exactly; to match its degrees of freedom, `civreghdfe` still counts an omitted column in K and in the residual and overidentification degrees of freedom when the fixed effects do not remove it exactly (it is absorbed only jointly by several unbalanced fixed effects, e.g. age = year - birth year with `absorb(id year)`), and absorbed exogenous regressors count among the instruments in the Cragg-Donald and Kleibergen-Paap Wald F degrees of freedom
- `partial()` variables are partialled out of the demeaned data with the estimation weights; collinear or absorbed `partial()` variables are dropped with a note. `e(partial)` and `e(partial_ct)` record them; without `absorb()`, `predict` returns only residuals (the full model's), as `ivreg2`
- Without `absorb()`, `predict, residuals` computes the residuals also outside `e(sample)`, as `ivreg2`'s `predict`; after `partial()` the partialling regression uses the estimation weights, where `ivreg2_p` ignores them
- Without `absorb()`, the constant is screened for collinearity first, so a regressor that is constant on the sample, or a full set of dummies, is omitted instead of `_cons`, as in `ivreg2`
- `coviv` gives LIML, Fuller and k-class the 2SLS covariance form (X'P_Z X)^-1 in every covariance estimator, as `ivreg2`

### Diagnostic-test options
- `orthog()`, `endogtest()`, `redundant()` and `partial()` accept factor-variable and time-series varlists, matched against the expanded model terms; repeated names count once and base levels are skipped
- A tested instrument or endogenous regressor that was dropped as collinear is an error (as in `ivreg2`)

### Degrees of Freedom
- Absorbed fixed effects consume degrees of freedom
- Apply the cluster-based degrees-of-freedom adjustment when standard errors are clustered

### Overidentification
- When instruments > endogenous variables, the model is overidentified
- Sargan/Hansen statistics are returned in `e(sargan)`, `e(sargan_df)`, and `e(sargan_p)` for supported overidentified specifications; see the installed help for inference details.

## Performance

The implementation uses the following methods:

**HDFE Absorption (shared with `creghdfe`):**
- **CG solver with symmetric Kaczmarz transform**: Forward and backward sweeps define the projection operator
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

Runtime depends on data and hardware; report the dataset, options, hardware, thread count, and software versions with timing comparisons.

## Comparison with ivreghdfe

| Feature | `ivreghdfe` | `civreghdfe` |
|---------|-------------|--------------|
| Implementation | Stata + Mata | C with OpenMP |
| HDFE Algorithm | Alternating projections | CG with Kaczmarz |
| Parallelization | Limited | Yes |

## See Also

- [ctools Overview](../../README.md)
- [creghdfe](creghdfe.md) - OLS regression with HDFE
- ivreghdfe - Original implementation (if installed)
- [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html) - Alternative IV estimation

## Validation and failure behavior

maxiter() and tolerance() must be positive. A fixed-effect projection that exhausts its iteration limit returns error 430. If memory allocation for the covariance matrix fails, the command returns an error.

# cppmlhdfe

`cppmlhdfe` estimates Poisson pseudo-maximum likelihood in C. It does not
require `ppmlhdfe` or `reghdfe`.
Models can include regressors, only fixed effects, or no absorbed effects.

```stata
cppmlhdfe y x1 x2, absorb(firm year) vce(cluster firm year) d(fe_sum)
predict double fitted, mu
cppmlhdfe y x, absorb(firm##c.trend region#year) exposure(population)
cppmlhdfe y x, vce(robust)
```

Numeric and string fixed effects and clusters, categorical interactions,
heterogeneous slopes (`g#c.z`, `g##c.z`, `g##c.(z t)`, and `g##(c.z c.t)` as in ppmlhdfe), factor regressors,
frequency/probability/analytic weights, and offset/exposure are supported.
The current implementation permits ten expanded FE terms and ten cluster
dimensions. Each cluster dimension must retain at least two clusters.

## Separation and estimation

The default `separation(fe simplex relu)` screens all-zero FE groups and
singleton chains, runs a native homogeneous simplex check, then uses an
iterated-rectifier (ReLU) check for separation from regressors and fixed effects,
including heterogeneous slopes. The separation checks are implemented in C. ReLU
uses weighted QR and follows `ppmlhdfe`'s projection and stopping rules.
`separation(mu)` enables `ppmlhdfe`'s fitted-mean screening and refitting rule;
`separation(all)` includes all four methods. `separation(none)` disables them. After removing separated observations, the
estimator rebuilds the sample, removes resulting singletons, and checks
collinearity again. The original regressor standardization is retained across
this refit, which matters for multiway covariance corrections.

IRLS preserves relative weights, including very small probability weights.
Response standardization uses centered variance. At tighter tolerances, convergence
also requires sufficiently small changes in the fitted log-means. Projection failure or exhausted IRLS
iterations returns 430 instead of posting a converged fit.

Fitted means are floored at 2^-46 (in units of sd(y)) for the IRLS weights
only, as in `ppmlhdfe`, and never capped. An IRLS step that would overflow a
fitted mean is halved until every mean is finite, as in `glm`'s IRLS; a
non-finite linear predictor stops the fit with error 430. When X'WX is
numerically singular far from the optimum, the step comes from a rank-revealing
QR of the partialled data. While some positive outcomes have fitted means
below the floor, the convergence test uses the deviance at the refitted linear
predictor (below) and also requires the fitted log-means to have settled to
the square root of `tolerance()`. The log likelihood and deviance are evaluated
from the linear predictor of the posted coefficients, with the fixed effects
recovered by one extra projection, the unit deviance in its cancellation-free
form and the saturated part by Stirling's series. When offsets or exposures vary greatly, the commands can differ in log likelihood,
deviance, and occasionally coefficients because of rounding in calculations
involving very small fitted means.

## Advanced controls

| Purpose | Options and defaults |
| --- | --- |
| Initial fitted means | `guess(simple)`, `guess(ols)`, `guess(variable varname)` |
| Scaling and collinearity | `standardize_data(1)`, `remove_collinear_variables(1)` |
| IRLS solves | `use_exact_solver(0)`, `use_exact_partial(1)`, `min_ok(1)` |
| Projection tolerance | `use_heuristic_tol(1)`, `start_inner_tol(1e-4)`, `itolerance(#)` |
| Step halving | `use_step_halving(0)`, `step_halving_memory(.9)`, `max_step_halving(2)` |
| ReLU | `relu_tol(1e-4)`, `relu_zero_tol(1e-8)`, `relu_maxiter(100)`, `relu_strict(0)`, `relu_accelerate(0)` |
| Simplex and mean screening | `simplex_tol(1e-12)`, `simplex_maxiter(1000)`, `mu_tol(1e-6)` |
| Projection solver | `accel(cg)`, `none`, `sd`, `aitken`, `hybrid`, or `lsmr` (any case, any leading part, `off` for `none`, as reghdfe reads them) |
| Projection transform | `transform(symmetric_kaczmarz)`, `kaczmarz`, or `cimmino` (any case, any leading part) |
| ppmlhdfe's GLM properties | `target_inner_tol(#)` (= `itolerance()`), `relu_debug(0)`, `relu_report_r2(0)` (= `r2`); all the properties above may also be given inside `absorb(..., )` by their full names |
| Column pools and LSMR | `poolsize(0)` (all columns), `btol(1e-8)`, `conlim(1e8)` |

`use_exact_partial(1)` (the default) projects the working response and the
original regressors each iteration. From the second iteration the projected
response is the working response less the previous iteration's fitted effects:
the residual is the same, but the projection tolerance applies to that
difference. Subtracting the previous iteration's fitted effects before projection reduces
projection error for zero outcomes whose fitted means approach zero. While a positive outcome has
its mean at the floor, whose working value is near y/floor, the working
response itself is projected. `use_exact_partial(0)` reuses the preceding
partialled working response and regressors. `use_exact_solver(1)` uses weighted pivoted QR. `min_ok()` requires
repeated qualifying convergence checks. Strict ReLU returns error 9010 if
it reaches the iteration limit without establishing separation. `tagsep(name)` switches
to a diagnostic-only ReLU run, overrides the requested separation methods, and
saves an indicator of separated observations without estimating a PPML model.
`zvarname(name)` saves the values used to establish separation. Rows excluded
from the diagnostic sample remain missing. These values may use a different
normalization from `ppmlhdfe`. With `r2`, as in `ppmlhdfe`, the certificate is
then regressed on the regressors and absorbed effects (with `creghdfe`, where
`ppmlhdfe` runs `reghdfe`) when observations were separated; its R-squared of 1
verifies it.

Individual effects can be saved with `absorb(firm_fe=firm year_fe=year)` or
`absorb(firm year, savefe)`. The latter creates `__hdfe1__`, `__hdfe2__`, etc.
For `trend_fe=firm##c.trend`, the slope is saved as `trend_feSlope1`; multiple
slopes use successive numbers. Their contributions reconstruct the saved
sum of effects. They are the IRLS-weighted least-squares fit of the saved sum on
the absorbed effects, solved by conjugate gradients from zero (minimum-norm
coefficients); `iterate()` bounds the iterations, and the stopping rule is
reghdfe's for saved effects, with the weighted root mean square of the sum in
place of its unit scale. If an automatically generated output name (`_ppmlhdfe_d` for `d2`,
`__hdfe1__` etc. for `savefe`) conflicts with a variable used by the model
(regressor, weight, offset, exposure, cluster or absorbed variable), the command
returns error 198 without changing the dataset. Individual coefficients are not uniquely identified and can
use a different normalization from the reference. Absorb suboptions also accept
`generate`, `keepsingletons`, `tolerance`, `iterate` (or reghdfe's
`maxiterations`), `dofadjustments`, `acceleration`, `transform`, `poolsize`,
`verbose(#)` (the fixed-effects object's verbosity: `verbose(-1)` hides the note
on dropped singletons), `prune`, `condition`, `rre()`, `algorithm()`,
`slopemethod()`, and the GLM properties by their full names, as `ppmlhdfe`
passes its remaining options to `reghdfe`. These partialling options leave the
estimates unchanged: `algorithm()` and `slopemethod()` are unused (as in
reghdfe 5's parser), `rre()` is accepted, `prune` is accepted without
`e(pruned)` (`ppmlhdfe`'s needs moremata), and `condition` (two absorbed effects
required) posts `e(finite_condition)` = -1, as `ppmlhdfe` does.
`precondition` and `residuals()` are refused, as by `ppmlhdfe` (which reports
Mata's code r(3498) for the first, where `cppmlhdfe` gives r(198), as for other
values reghdfe's parser rejects with r(3598)); `groupvar()`, which `ppmlhdfe`
accepts, is not supported.

## Output and replay

The output follows `ppmlhdfe`: the notes on dropped singletons and on
observations separated by a fixed effect (`(dropped N observations that are
either singletons or separated by a fixed effect)`), on the separated
observations each method drops (`(simplex method dropped ...)`, `(ReLU method
dropped ... in N iterations)`, `(IRLS step detected ...)`), the warnings for
`keepsingletons` and for an outcome with very low values after standardizing,
the IRLS log, `ppmlhdfe`'s header (observations, residual df, Wald test,
deviance, log pseudolikelihood, pseudo R2, `e(title3)` and the clusters per
dimension), the coefficient table, and reghdfe's table of absorbed degrees of
freedom (not without absorbed effects). `cppmlhdfe` without a variable list
replays the last fit with the same parts; `eform` (or its synonym `irr`, both
shown as exp(b)), `noheader`, `notable` and the display options apply, and the
replay returns `_coef_table`'s `r()`. `cppmlhdfe, version` shows the command's
location and version and the ctools plugin; `cppmlhdfe, reload` tells how to
reinstall ctools (`ctools, update`), where `ppmlhdfe`'s reinstalls from its
author's working copy.

## Covariance, inference, and prediction

Robust and clustered covariance include the reported constant. Multiway
clustering uses inclusion-exclusion over all nonempty cluster subsets, with a
multiplier based on the minimum marginal cluster count. Negative eigenvalues of
the full standardized covariance are set to zero before rescaling, as in
`ppmlhdfe`. Cost grows as `2^Q - 1`. Single-way calculations use a centered slope
system to avoid cancellation when regressor means are large.

`vce(bootstrap [, bootopts])` gives pairs-bootstrap standard errors, computed
in the plugin (ppmlhdfe has no such option), and posts what `vce(bootstrap)`
posts after official commands, so `estat bootstrap` works; the fit's robust
covariance is `e(V_modelbased)`. The options are creghdfe's: `reps(#)` (50),
`seed(#)`, `cluster(clustvar)`, `strata(varlist)`, `draws(stata|fast)` (default `stata`), `mse`,
`level(#)`, `nodots`/`dots(#)` and `saving()`; `e(p)` is the p-value of the Wald `e(chi2)`,
and `mse` and `level()` behave as in creghdfe. Each replicate refits PPML on the drawn rows
weighted by their multiplicities: screening on the draws, the separation checks
again (rows they find leave the replicate, which is refit), and IRLS started
from the full fit. With `draws(stata)` the results reproduce
`bootstrap _b, reps(#) seed(#): cppmlhdfe ...` to the IRLS tolerance. Weights
are allowed: aweights and pweights multiply the resampled rows; with
`cluster()` each drawn copy of an fweighted row keeps its weight; without it,
fweights are resampled as the observations they stand for, as in `creghdfe`
(with `draws(stata)`, the prefix's draws on the expanded data). Heterogeneous
slopes are re-estimated in each replicate. Not with `tagsep()` or without
regressors.

`fevce(bootstrap [, feopts])` gives bootstrap standard errors of the saved
fixed effects (`savefe`, `newvar=`) and of `d()`. Each replicate recovers its
effects as the fit does (minimum-norm coefficients of its d); a level's SE is
the standard deviation over the replicates in which it is estimable (drawn,
and its full-sample connected set over all effects not split). Names:
`newvar_se`, `__fese#__` for `__hdfe#__` (replaced on reruns), and `d_se` for
`d()` unless `dse()` names it; `sestub()`, `nreps` (`..._nrep`,
`__fenrep#__`) and `ci` (`..._lb`/`..._ub`, `__felb#__`/`__feub#__`) as in
creghdfe. A slope saved in `__hdfe#__Slope1` has its SE in `__fese#__Slope1`;
it is estimable where the replicate's slope variable varies (and where the
full fit estimates the slope but a replicate cannot, that replicate does not
estimate the level's intercept either). With `vce(bootstrap)` both share the
replicates, and `level()`, `mse` and `dots()` are those of `vce()`. With two or more absorbed
intercepts a normalization warning is shown and `e(fe_norm_warning)` is 1.
Results: `e(fevce)`, `e(fe_se)`, `e(fe_nrep)`, `e(fe_lb)`, `e(fe_ub)`,
`e(fe_dse)`, `e(fe_bs_table)`, `e(fe_N_reps)`, `e(fe_split_reps)`,
`e(fe_rngstate)`.

`wildbootstrap[(wildopts)]` adds score bootstrap p-values (Kline and Santos,
as boottest computes them for ML): per tested coefficient the model is refit
under H0: b = 0 (the regressor moves into the offset), and each replicate
perturbs the clusters' efficient score contributions and studentizes them,
recentered by cluster shares. e(b), e(V) and the analytic table stay. It needs
`vce(robust)` or one-way `vce(cluster)`. Options as in creghdfe, plus `ci` for
intervals by test inversion (a refit per trial value, about 30-50 per
coefficient); `nonull` is not available. Unweighted, the p-values equal
boottest's after `poisson` with dummies for the same seed. With weights the
score contributions are weighted (boottest leaves them unweighted, and its
fweighted results then differ from its results on the expanded data, which
cppmlhdfe reproduces). Combines with `fevce(bootstrap)`, not with
`vce(bootstrap)`; `nowildtable` hides the table. Results: `e(wb_method)` =
`score`, `e(wb_t)` (Wald z), `e(wb_score)` (the score statistic the bootstrap
studentizes, shown in the table), `e(wb_p)`, `e(wb_fail)`, `e(wb_ci)` and `e(wb_cistatus)`
with `ci`, `e(wb_reps)`, `e(wb_misreps)`, `e(wb_enumerated)`, `e(wb_N_clust)`,
`e(wb_rngstate)`.

Inference uses normal critical values and a Wald chi-squared test. Residual
degrees of freedom are stored as `e(df)`. The command also reports absorbed-DF
components, nesting, singleton/separation counts, and `e(dof_table)`.
`e(title)` is `HDFE PPML regression`, or `PPML regression` without absorbed
effects (plain text, without `ppmlhdfe`'s help link), and `e(marginsnotok)` is
`ppmlhdfe`'s list. `e(rss)` and `e(rmse)` are `ppmlhdfe`'s statistics of the
final working regression (squared partialled working residuals weighted by the
IRLS weights rescaled to sum to N, or to the fweight total; the root mean square
on N - df_a - rank - df_a_nested). `e(ic2)` counts the partialling
subiterations of the IRLS; it equals `ppmlhdfe`'s with one absorbed effect or
none, and differs with several because the projection algorithms differ.
`e(title2)` and `e(N_hdfe)` count the `absorb()` terms, as `ppmlhdfe`
(`absorb(id##c.x year)` is two); `e(title3)` is `Statistics robust to
heteroskedasticity` with clustering; `e(indepvars)` lists the estimated
regressors as `ppmlhdfe` does (factor levels in `#bn.` form, without base levels
and omitted terms) followed by `_cons`. `e(footnote)` and `e(estat_cmd)` name
`_cppmlhdfe_footnote` and `_cppmlhdfe_estat` (`ppmlhdfe` names `reghdfe`'s
programs, which refuse other commands); `estat` offers the same subcommands,
`summarize`, `vce` and `ic` with `e(df_m) + 1` degrees of freedom. With a
constant outcome, the omitted
regressors keep their names in `e(b)` where `ppmlhdfe` prefixes them with `o.`;
`e(V)` and the output agree.

`d(name)` saves the sum of fixed effects; bare `d` uses `_ppmlhdfe_d`.
`predict` supports `xb`, `mu` (default, with `ppmlhdfe_p`'s note "(option mu
assumed; predicted mean of depvar)"), `eta`/`xbd`, `d`, `stdp`, `response`,
`scores`, `pearson`, `deviance`, `anscombe`, and `working`, including after
`estimates store`/`estimates restore`. As in the reference, absorbed models need
`d()` for predictions other than `xb`. `xb` includes the constant and offset or
log exposure, and excludes the absorbed effects. The new variables carry
`ppmlhdfe_p`'s labels (`Predicted mean of y`, `Linear prediction: xb + d[...]`,
`Pearson residual`, ...). Without absorbed effects, `d` gives 0 (`ppmlhdfe_p`
fails there).

## Validation and compatibility limits

`validation/stata/cppmlhdfe.do` includes the original examples,
the named cases `cppmlhdfe_precision_cases`, `cppmlhdfe_compatibility_cases`,
`cppmlhdfe_upstream_cases`, and `cppmlhdfe_advanced_cases` in that same file. These compare
actual retained rows, coefficient
names, full covariance including the constant, standard errors, inference,
absorbed degrees of freedom, and predictions. Coverage includes randomized
weights/clusters/offsets, disconnected FE graphs, separation, string groups,
response/weight scaling, heterogeneous slopes, and upstream difficult examples.

Default fits require relative SE error below `1e-5`. Precision comparisons use
`tolerance(1e-14)` and exact partialling/solves in the reference, with six
significant digits required. Prediction comparisons use three converged reference
iterations because its saved FE predictions can otherwise lag the coefficient
solution. Native separation/projection tests run under ASan and UBSan.

Run `python3 validation/run.py cppmlhdfe` for the current Stata checks.
The suite covers randomized default fits, comparisons of coefficients, standard errors, and covariance matrices
using the tolerances stated above, weighted heterogeneous slopes, saved-effect reconstruction,
diagnostic-only ReLU, initialization, and alternative solvers. Native checks run
through `python3 validation/run.py --native`; set
`CTOOLS_SANITIZERS=address,undefined` to check memory and undefined behavior.

The upstream `slopes4` regression check requires the reference's default
sample: 14 retained observations, three separated rows, and ten singletons.
The implementation uses joint intercept/slope projection, warm residual
reuse, and tolerance-dependent stopping rules. For a saturated upstream model,
the reference is initialized at its exact mean to avoid premature convergence
from rounded zero deviance. Multiple absorbed slopes are checked against the
equivalent explicit-regressor reference model because the installed release's
absorbed-slope preconditioner errors on that fixture. Alternative native solvers
are also checked against the reference CG solution; the installed reference's
hybrid solver fails internally on the separated test sample.

The tests compare final results for the listed models. Intermediate calculations
may differ. Unsupported `ppmlhdfe` options are listed below.
Mobility-group output (`groupvar`) and hidden acceleration scheduling controls
are not implemented. Unsupported options produce an explicit error. Legacy
`septolerance()` is accepted but does not set the `mu` screening threshold;
use `mu_tol()` when explicitly enabling that method.

Run the Stata tests separately on each supported platform. CI builds all four
plugins and runs C tests, but does not run Stata. Keep ado files and plugins together; the command
rejects an older plugin API.

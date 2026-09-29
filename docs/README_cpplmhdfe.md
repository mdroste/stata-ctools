# cpplmhdfe

`cpplmhdfe` estimates Poisson pseudo-maximum likelihood in C, with no runtime
fallback to `ppmlhdfe` or `reghdfe`. The public command is spelled `cpplmhdfe`.
Models can include regressors, only fixed effects, or no absorbed effects.

```stata
cpplmhdfe y x1 x2, absorb(firm year) vce(cluster firm year) d(fe_sum)
predict double fitted, mu
cpplmhdfe y x, absorb(firm##c.trend region#year) exposure(population)
cpplmhdfe y x, vce(robust)
```

Numeric and string fixed effects and clusters, categorical interactions,
heterogeneous slopes (`g#c.z`, `g##c.z`, `g##c.(z t)`), factor regressors,
frequency/probability/analytic weights, and offset/exposure are supported.
The current implementation permits ten expanded FE terms and ten cluster
dimensions. Each cluster dimension must retain at least two clusters.

## Separation and estimation

The default `separation(fe simplex relu)` screens all-zero FE groups and
singleton chains, runs a native homogeneous simplex check, then uses an
iterated-rectifier (ReLU) check for separation from regressors and fixed effects,
including heterogeneous slopes. These are separate C implementations. ReLU
uses weighted QR and reference-compatible projection and stopping rules.
`separation(mu)` enables the reference's fitted-mean screening and refit rule;
`separation(all)` includes all four methods. `separation(none)` disables them. After removing separated observations, the
estimator rebuilds the sample, removes resulting singletons, and checks
collinearity again. The original regressor standardization is retained across
this refit, which matters for multiway covariance corrections.

IRLS preserves relative weights, including very small probability weights.
Response standardization uses centered variance. Tight fits additionally check
stability of identified fitted log-means. Projection failure or exhausted IRLS
iterations returns 430 instead of posting a converged fit.

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
| Projection solver | `accel(cg)`, `none`, `sd`, `aitken`, `hybrid`, or `lsmr` |
| Projection transform | `transform(symmetric_kaczmarz)`, `kaczmarz`, or `cimmino` |
| Column pools and LSMR | `poolsize(0)` (all columns), `btol(1e-8)`, `conlim(1e8)` |

`use_exact_partial(0)` reuses the preceding partialled working response and
regressors. `use_exact_solver(1)` uses weighted pivoted QR. `min_ok()` requires
repeated qualifying convergence checks. Strict ReLU returns error 9010 when
its iteration limit is reached without a certificate. `tagsep(name)` switches
to a diagnostic-only ReLU run, overrides the requested separation methods, and
saves the separation indicator without estimating a PPML model. Add
`zvarname(name)` to save the certificate. Rows excluded from the diagnostic
sample remain missing. Certificates need not use the same arbitrary scale as
the reference.

Individual effects can be saved with `absorb(firm_fe=firm year_fe=year)` or
`absorb(firm year, savefe)`. The latter creates `__hdfe1__`, `__hdfe2__`, etc.
For `trend_fe=firm##c.trend`, the slope is saved as `trend_feSlope1`; multiple
slopes use successive numbers. Their contributions reconstruct the saved
sum of effects. Individual coefficients are not uniquely identified and can
use a different normalization from the reference. Absorb suboptions also accept
`generate`, `keepsingletons`, `tolerance`, `iterate`, `dofadjustments`,
`acceleration`, `transform`, and `poolsize`.

## Covariance, inference, and prediction

Robust and clustered covariance include the reported constant. Multiway
clustering uses inclusion-exclusion over all nonempty cluster subsets, with a
multiplier based on the minimum marginal cluster count. Negative eigenvalues of
the full standardized covariance are set to zero before rescaling, as in
`ppmlhdfe`. Cost grows as `2^Q - 1`. Single-way calculations use a centered slope
system to avoid cancellation when regressor means are large.

Inference uses normal critical values and a Wald chi-squared test. Residual
degrees of freedom are stored as `e(df)`. The command also reports absorbed-DF
components, nesting, singleton/separation counts, and `e(dof_table)`.

`d(name)` saves the sum of fixed effects; bare `d` uses `_ppmlhdfe_d`.
`predict` supports `xb`, `mu` (default), `eta`/`xbd`, `d`, `stdp`, `response`,
`scores`, `pearson`, `deviance`, `anscombe`, and `working`, including after
`estimates store`/`estimates restore`. As in the reference, absorbed models need
`d()` for predictions other than `xb`. `xb` includes the constant and offset or
log exposure, and excludes the absorbed effects.

## Validation and compatibility limits

`validation/validate_cpplmhdfe.do` includes the original examples,
`cpplmhdfe_precision_cases.do`, `cpplmhdfe_compatibility_cases.do`,
`cpplmhdfe_upstream_cases.do`, and `cpplmhdfe_advanced_cases.do`. These compare
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

The final macOS ARM run passed **1,065 PPML checks** and **23 shared regression
checks**, with no skips, against installed `ppmlhdfe` 2.3.3 (02nov2025).
Across 360 randomized default fits, the largest relative slope-SE difference
was `4.630861419e-7`. All 450 stringent randomized coefficient/SE/covariance
comparisons passed. The 40 advanced-option checks include weighted
heterogeneous slopes, saved-effect reconstruction, diagnostic-only ReLU,
initialization, and alternative solvers. Native ASan/UBSan tests and the plugin
dependency-contract check also passed.

The pathological upstream `slopes4` fixture now matches the reference's default
sample: 14 retained observations, three separated rows, and ten singletons.
The correction reproduces its joint intercept/slope projection, warm residual
reuse, and tolerance-dependent stopping rule. For a saturated upstream model,
the reference is initialized at its exact mean to avoid premature convergence
from rounded zero deviance. Multiple absorbed slopes are checked against the
equivalent explicit-regressor reference model because the installed release's
absorbed-slope preconditioner errors on that fixture. Alternative native solvers
are also checked against the reference CG solution; the installed reference's
hybrid solver fails internally on the separated test sample.

This establishes numerical compatibility for the tested models, not identical
floating-point iteration paths or support for every hidden reference option.
Reference debugging/benchmarking options (`r2`, `rre`), mobility-group output
(`groupvar`), graph pruning (`prune`), and hidden acceleration scheduling controls
are not implemented. Unsupported options produce an explicit error. Legacy
`septolerance()` is accepted but does not set the `mu` screening threshold;
use `mu_tol()` when explicitly enabling that method.

Validation and the rebuilt plugin in this checkout use macOS ARM. Other platforms
must be rebuilt and tested separately. Keep the ado files and plugin together;
the command rejects an older plugin API.

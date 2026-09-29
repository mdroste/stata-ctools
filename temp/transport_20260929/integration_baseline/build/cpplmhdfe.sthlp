{smcl}
{* *! version 1.1.0 25Sep2026}{...}
{viewerjumpto "Syntax" "cpplmhdfe##syntax"}{...}
{viewerjumpto "Description" "cpplmhdfe##description"}{...}
{viewerjumpto "Options" "cpplmhdfe##options"}{...}
{viewerjumpto "Examples" "cpplmhdfe##examples"}{...}
{viewerjumpto "Stored results" "cpplmhdfe##results"}{...}
{title:Title}

{phang}
{bf:cpplmhdfe} {hline 2} C-accelerated Poisson pseudo-maximum likelihood regression with high-dimensional fixed effects


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cpplmhdfe}
{depvar}
[{indepvars}]
{ifin}
[{it:weight}]
{cmd:,}
[{opt a:bsorb(absvars)} {it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(varlist)}}categorical interactions and heterogeneous slopes; optional{p_end}
{synopt:{opt exp:osure(varname)}}exposure variable (enters as log offset){p_end}
{synopt:{opt off:set(varname)}}offset variable (enters linearly in the linear predictor){p_end}
{synopt:{opt d(newvar)}}save the sum of fixed effects for prediction{p_end}
{synopt:{opt keepsin:gletons}}retain singleton observations{p_end}
{synopt:{opt dof:adjustments(doftype)}}degrees of freedom adjustment method{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt robust} or {opt cluster} {it:clustvars}{p_end}

{syntab:IRLS Convergence}
{synopt:{opt irlstol:erance(#)}}IRLS convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt irlsmax:iter(#)}}maximum IRLS iterations; default is {cmd:10000}{p_end}
{synopt:{opt septol:erance(#)}}legacy option; small means do not establish separation{p_end}

{syntab:CG Solver}
{synopt:{opt tol:erance(#)}}deviance convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt iter:ate(#)}}maximum CG iterations; default is {cmd:10000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Reporting}
{synopt:{opt v:erbose(#)}}display detailed timing and convergence information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{cmd:aweight}s, {cmd:fweight}s, and {cmd:pweight}s are allowed; see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cpplmhdfe} estimates Poisson pseudo-maximum likelihood (PPML) regressions
with multiple high-dimensional fixed effects using an IRLS (Iteratively Reweighted
Least Squares) algorithm. It is a C-accelerated replacement for {cmd:ppmlhdfe}.

{pstd}
PPML is the standard estimator for gravity models in trade economics and for count
data with many fixed effects. The estimator is consistent for the conditional mean
E[y|x] = exp(x'b) without requiring the data to follow a Poisson distribution
(Santos Silva and Tenreyro, 2006).

{pstd}
The dependent variable must be non-negative but need not be an integer.

{pstd}
The algorithm uses the ctools HDFE infrastructure (conjugate gradient solver,
singleton detection) inside an IRLS outer loop. At each iteration, Poisson
working weights (mu * user_weights) are used for weighted FE projection via
the CG solver.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt absorb(varlist)} specifies one or more categorical variables whose effects
are to be absorbed (partialled out). Numeric and string categories are allowed.
Examples include {cmd:g#h}, {cmd:g#c.z}, {cmd:g##c.z}, and
{cmd:g##c.(z t)}. Up to ten expanded FE terms are supported. Omitting
{cmd:absorb()} estimates a model with a constant and no absorbed effects.
Save individual effects with {cmd:absorb(firm_fe=firm year_fe=year)} or
{cmd:absorb(firm year, savefe)}. The latter uses names {cmd:__hdfe1__}, etc.
Slope coefficients use suffixes {cmd:Slope1}, {cmd:Slope2}, etc.
Their contributions reconstruct the sum saved by {cmd:d(newvar)}; individual
normalizations need not equal the reference. Absorb suboptions also accept
{cmd:generate}, {cmd:keepsingletons}, {cmd:tolerance()}, {cmd:iterate()},
{cmd:dofadjustments()}, {cmd:acceleration()}, {cmd:transform()}, and
{cmd:poolsize()}. Models without regressors are allowed.

{phang}
{opt exposure(varname)} specifies an exposure variable. This adds the logarithm of {it:varname} as the offset and is commonly used in rate models.

{phang}
{opt offset(varname)} specifies an offset variable that enters the linear
predictor additively: eta = X*beta + FE + offset.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies the type of standard error. Options are
{opt robust} for Eicker-Huber-White sandwich standard errors and
{opt cluster} {it:clustvars} for one-way or multiway cluster-robust standard errors
(up to ten numeric or string cluster variables). Multiway covariance uses
inclusion-exclusion over cluster intersections and the minimum marginal cluster
count for its finite-sample multiplier. Negative eigenvalues of the full
covariance are set to zero, matching {cmd:ppmlhdfe}. Every cluster dimension must
retain at least two clusters. Computation grows with the number of nonempty
cluster subsets (2^Q-1).

{pstd}
Inference uses normal critical values and Wald chi-squared tests, as in
{cmd:ppmlhdfe}. The normalized constant and its covariance are included in {cmd:e(b)} and {cmd:e(V)}.

{dlgtab:IRLS Convergence}

{phang}
{opt irlstolerance(#)} specifies the convergence criterion for the IRLS
loop. Convergence is declared when the relative change in deviance falls
below this threshold: |dev - dev_old| / max(min(dev, dev_old), 0.1) < tol,
also requiring sufficient accuracy in the FE projection. Explicit tolerances
below 1e-10 additionally require changes in fitted log-means below 1e-8. Deviance is evaluated
on the standardized response with a common weight normalization. Default is
{cmd:tolerance()}, or 1e-8. For demanding numerical comparisons, specify
{cmd:tolerance(1e-14)} in both commands; their default iteration paths can differ.

{phang}
{opt irlsmaxiter(#)} specifies the maximum number of IRLS iterations.
Default is 10000. {cmd:maxiterations()} is an alias.

{phang}
{opt separation(fe simplex relu)} is the default. All-zero FE groups and
singletons are screened first. Native simplex and iterated-rectifier (ReLU)
checks then detect general regressor and FE separation. The model is rebuilt
after deleting separated rows, including another singleton and collinearity
check. {cmd:separation(mu)} enables the reference fitted-mean screening rule;
{cmd:separation(all)} runs all four methods. {cmd:separation(none)} disables checks.

{phang}
{opt septolerance(#)} is a legacy compatibility option. Small fitted means alone
do not activate mean screening by themselves. Explicit {cmd:separation(mu)}
uses {cmd:mu_tol()} instead. All separation methods run in C without calling
{cmd:ppmlhdfe}.

{phang}
{opt d(newvar)} saves the normalized sum of fixed effects for prediction.
Bare {cmd:d} uses {cmd:_ppmlhdfe_d}. {cmd:keepsingletons} disables singleton
removal.

{dlgtab:Advanced controls}

{phang}
{cmd:guess(simple)} (default), {cmd:guess(ols)}, and
{cmd:guess(variable varname)} select initial fitted means.
{cmd:standardize_data(1)} and {cmd:remove_collinear_variables(1)} control
standardization and the preliminary collinearity check.

{phang}
{cmd:use_exact_solver(1)} selects weighted pivoted QR; default is 0.
{cmd:use_exact_partial(1)} (default) projects the original working response and
regressors each iteration; 0 reuses preceding residuals.
{cmd:min_ok(1)} sets the required number of qualifying convergence checks.
{cmd:use_heuristic_tol(1)}, {cmd:start_inner_tol(1e-4)}, and {cmd:itolerance()}
control adaptive projection accuracy.

{phang}
{cmd:use_step_halving(0)}, {cmd:step_halving_memory(.9)}, and
{cmd:max_step_halving(2)} control damping after deviance increases.

{phang}
ReLU controls are {cmd:relu_tol(1e-4)}, {cmd:relu_zero_tol(1e-8)},
{cmd:relu_maxiter(100)}, {cmd:relu_strict(0)}, and {cmd:relu_accelerate(0)}.
Strict mode returns 9010 when the ReLU limit is exhausted.
{cmd:simplex_tol(1e-12)}, {cmd:simplex_maxiter(1000)}, and {cmd:mu_tol(1e-6)}
control the other separation methods. {cmd:tagsep(newvar)} switches to a
ReLU diagnostic-only run, overriding the requested separation methods and
saving the indicator without estimating a model. Add {cmd:zvarname(newvar)}
to save the certificate. Excluded observations remain missing; certificate
scales can differ from the reference.

{phang}
{cmd:acceleration()} accepts {cmd:cg} (default), {cmd:none}, {cmd:sd},
{cmd:aitken}, {cmd:hybrid}, and {cmd:lsmr}. {cmd:transform()} accepts
{cmd:symmetric_kaczmarz} (default), {cmd:kaczmarz}, and {cmd:cimmino}.
{cmd:poolsize(0)} projects all columns together; positive values limit each pool.
LSMR accepts {cmd:btol(1e-8)} and {cmd:conlim(1e8)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse ships}{p_end}

{pstd}Basic PPML with fixed effects{p_end}
{phang2}{cmd:. cpplmhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, absorb(ship) vce(robust)}{p_end}

{pstd}With exposure variable{p_end}
{phang2}{cmd:. cpplmhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, absorb(ship) exposure(service) vce(robust)}{p_end}

{pstd}Two-way fixed effects{p_end}
{phang2}{cmd:. cpplmhdfe accident op_75_79, absorb(ship type) vce(cluster ship)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cpplmhdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df)}}residual degrees of freedom (not used for t inference){p_end}
{synopt:{cmd:e(ll)}}log pseudolikelihood{p_end}
{synopt:{cmd:e(ll_0)}}log pseudolikelihood of null model{p_end}
{synopt:{cmd:e(deviance)}}Poisson deviance{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared (McFadden){p_end}
{synopt:{cmd:e(chi2)}}Wald chi-squared statistic{p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster dimensions{p_end}
{synopt:{cmd:e(N_clust1)}, etc.}marginal cluster counts{p_end}
{synopt:{cmd:e(ic)}}number of IRLS iterations{p_end}
{synopt:{cmd:e(converged)}}1 if IRLS converged, 0 otherwise{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(num_singletons)}}number of singleton observations dropped{p_end}
{synopt:{cmd:e(num_separated)}}number of separated observations dropped{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom absorbed by FEs{p_end}
{synopt:{cmd:e(N_clust)}}minimum marginal cluster count (if clustering){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cpplmhdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(absorb)}}absorbed variables{p_end}
{synopt:{cmd:e(vce)}}VCE type{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustering){p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Santos Silva, J.M.C. and Tenreyro, S. (2006).
"The Log of Gravity."
{it:Review of Economics and Statistics} 88(4): 641-658.

{phang}
Correia, S., Guimaraes, P., and Zylkin, T. (2020).
"Fast Poisson estimation with high-dimensional fixed effects."
{it:Stata Journal} 20(1): 95-115.
{p_end}


{marker author}{...}
{title:Author}

{pstd}
Part of the {browse "https://github.com/mdroste/stata-ctools":ctools} package.
{p_end}

{title:Prediction and estimation sample}

{pstd}
{cmd:predict newvar, xb} returns Xb including the constant and the offset or
log exposure. The default {cmd:mu} returns the fitted mean. Other supported
options are {cmd:eta}, {cmd:xbd}, {cmd:d}, {cmd:stdp}, {cmd:response}, {cmd:scores},
{cmd:pearson}, {cmd:deviance}, {cmd:anscombe}, and {cmd:working}. For absorbed
models, options other than {cmd:xb} require {cmd:d()} at estimation. These
predictions also work after {cmd:estimates store}/{cmd:estimates restore}.

{pstd}
{cmd:e(sample)} marks the retained observations after missing values, singleton
chains, and FE/regressor separation are removed. With fweights, {cmd:e(N)} is the
sum of frequency weights rather than the number of marked physical rows.
Weight expressions such as {cmd:[aw=2*w]} are evaluated into temporary doubles.

{title:Validation and failure behavior}

{pstd}
Fixed-effect and IRLS iteration limits and tolerances must be positive. A nonconverged projection or exhausted IRLS fit returns error 430. Joint Wald statistics use the retained regressor indices, so moving omitted columns does not change the reported test.
{p_end}

{pstd}
Numerical tests compare coefficients, full covariance, retained samples, and
predictions with {cmd:ppmlhdfe}. The pathological upstream heterogeneous-slope
fixture now matches its default sample: 14 retained rows, three separated rows,
and ten singletons. Iteration paths and hidden debugging options are not all
identical. See {cmd:docs/README_cpplmhdfe.md} for validation controls and limits.

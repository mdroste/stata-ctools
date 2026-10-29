{smcl}
{* *! version 1.2.0 20261029}{...}
{viewerjumpto "Syntax" "cppmlhdfe##syntax"}{...}
{viewerjumpto "Description" "cppmlhdfe##description"}{...}
{viewerjumpto "Options" "cppmlhdfe##options"}{...}
{viewerjumpto "Prediction" "cppmlhdfe##prediction"}{...}
{viewerjumpto "Examples" "cppmlhdfe##examples"}{...}
{viewerjumpto "Stored results" "cppmlhdfe##results"}{...}
{title:Title}

{phang}
{bf:cppmlhdfe} {hline 2} C-accelerated Poisson pseudo-maximum likelihood regression with high-dimensional fixed effects


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cppmlhdfe}
{depvar}
[{indepvars}]
{ifin}
[{it:weight}]
[{cmd:,} {opt a:bsorb(absvars)} {it:options}]

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
{synopt:{opt sep:aration(methods)}}separation checks; default {cmd:fe simplex relu}{p_end}
{synopt:{opt septol:erance(#)}}compatibility option; has no effect{p_end}

{syntab:Projection and convergence}
{synopt:{opt tol:erance(#)}}deviance convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt iterate(#)}}maximum CG iterations; default is {cmd:10000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Reporting}
{synopt:{opt v:erbose(#)}}display detailed timing and convergence information{p_end}
{synopt:{opt ef:orm}}report exponentiated coefficients{p_end}
{synopt:{opt ir:r}}report incidence-rate ratios; {opt eform} with the column titled IRR{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{cmd:aweight}s, {cmd:fweight}s, and {cmd:pweight}s are allowed; see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cppmlhdfe} estimates Poisson pseudo-maximum likelihood (PPML) regressions
with multiple high-dimensional fixed effects using an IRLS (Iteratively Reweighted
Least Squares) algorithm. The fixed-effect projections and likelihood calculations run in C.

{pstd}
The dependent variable must be nonnegative but need not be an integer. The
conditional mean is exp(Xb + fixed effects + offset). Factor-variable and
time-series notation are allowed. {cmd:e(sample)} identifies the retained
observations after missing values, singleton groups, and separation are removed.


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
Their contributions reconstruct the sum saved by {cmd:d(newvar)}; individual fixed-effect levels depend on normalization. Absorb suboptions also accept
{cmd:generate}, {cmd:keepsingletons}, {cmd:tolerance()}, {cmd:iterate()},
{cmd:dofadjustments()}, {cmd:acceleration()}, {cmd:transform()}, and
{cmd:poolsize()}. Models without regressors are allowed.

{phang}
{opt exposure(varname)} specifies an exposure variable. This adds the logarithm of {it:varname} to the linear predictor. Values must
be positive. It cannot be combined with {opt offset()}.

{phang}
{opt offset(varname)} specifies an offset variable that enters the linear
predictor additively: eta = X*beta + FE + offset.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies the type of standard error. The default is {opt robust}. Options are {opt robust} for Eicker-Huber-White sandwich standard errors and
{opt cluster} {it:clustvars} for one-way or multiway cluster-robust standard errors
(up to ten numeric or string cluster variables). Multiway covariance uses
inclusion-exclusion over cluster intersections and the minimum marginal cluster
count for its finite-sample multiplier. Negative eigenvalues of the full
covariance are set to zero, matching {cmd:ppmlhdfe}. Every cluster dimension must
retain at least two clusters. Computation grows with the number of nonempty
cluster subsets (2^Q-1). {cmd:vce(unadjusted)} is accepted as an alias for
{cmd:vce(robust)}; it does not select a model-based Poisson covariance.

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
{cmd:tolerance()}, or 1e-8.

{phang}
{opt irlsmaxiter(#)} specifies the maximum number of IRLS iterations.
Default is 10000. {cmd:maxiterations()} is an alias.

{phang}
{opt separation(fe simplex relu)} is the default. All-zero FE groups and
singletons are screened first. Native simplex and iterated-rectifier (ReLU)
checks then detect general regressor and FE separation. The model is rebuilt
after deleting separated rows, including another singleton and collinearity
check. {cmd:separation(mu)} enables fitted-mean screening;
{cmd:separation(all)} runs all four methods. {cmd:separation(none)} disables checks.

{phang}
{opt septolerance(#)} is accepted for compatibility and has no effect on
separation. Use {opt mu_tol()} with {cmd:separation(mu)} to control
fitted-mean screening.


{phang}
{opt d(newvar)} saves the normalized sum of fixed effects for prediction.
{opt d2} uses {cmd:_ppmlhdfe_d}, replacing that variable if it exists. {cmd:keepsingletons} disables singleton
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
to save the certificate. Excluded observations remain missing.

{phang}
{cmd:acceleration()} accepts {cmd:cg} (default), {cmd:none}, {cmd:sd},
{cmd:aitken}, {cmd:hybrid}, and {cmd:lsmr}. {cmd:transform()} accepts
{cmd:symmetric_kaczmarz} (default), {cmd:kaczmarz}, and {cmd:cimmino}.
{cmd:poolsize(0)} projects all columns together; positive values limit each pool.
LSMR accepts {cmd:btol(1e-8)} and {cmd:conlim(1e8)}.


{phang}
{opt tolerance(#)} sets the deviance convergence tolerance; default 1e-8.
{opt irlstolerance()} overrides it for the outer iteration. {opt iterate(#)}
sets the maximum number of inner projection iterations; default 10000.
Iteration limits and tolerances must be positive. An exhausted projection or
IRLS loop returns error 430.

{phang}
{opt dofadjustments(list)} controls absorbed degrees of freedom. The default
{cmd:all} applies pairwise redundancy checks, cluster nesting adjustments, and
checks for redundant slope terms. Other choices are {cmd:none}, or a combination
of {cmd:firstpair} or {cmd:pairwise}, {cmd:clusters}, and {cmd:continuous}.

{phang}
{opt noabsorb} explicitly requests a model without absorbed effects and cannot
be combined with {opt absorb()}. {opt cluster(varlist)} is an alias for
{cmd:vce(cluster varlist)}.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.
{opt verbose(#)} sets reporting detail; default 0. {opt timeit} requests timing
output. {opt nolog} is accepted; the command does not print an iteration log.

{phang}
{opt eform} reports exponentiated coefficients. {opt irr} gives the same values
with the column labeled IRR. Stored coefficients remain on the log scale.


{marker prediction}{...}
{title:Prediction and estimation sample}

{pstd}
{cmd:predict newvar, xb} returns Xb including the constant and the offset or
log exposure. The default {cmd:mu} returns the fitted mean. Other supported
options are {cmd:eta}, {cmd:xbd}, {cmd:d}, {cmd:stdp}, {cmd:response}, {cmd:scores},
{cmd:pearson}, {cmd:deviance}, {cmd:anscombe}, and {cmd:working}. For absorbed
models, options other than {cmd:xb} require {opt d()} or {opt d2} at estimation.
Keep the saved fixed-effect, offset, and exposure variables for predictions
after {cmd:estimates store}/{cmd:estimates restore}.

{pstd}
Let mu denote the fitted mean and eta its logarithm. {opt eta} and {opt xbd}
return eta, including absorbed effects; {opt d} returns the saved sum of fixed
effects. {opt stdp} is the standard error of the linear predictor excluding
absorbed-effect uncertainty. {opt response} and {opt scores} return y-mu;
{opt pearson} returns (y-mu)/sqrt(mu); {opt working} returns (y-mu)/mu.
{opt deviance} returns 2*(mu-y+y*ln(y/mu)), with the log term set to zero
when y=0. {opt anscombe} returns 1.5*(y^(2/3)-mu^(2/3))/mu^(1/6).

{pstd}
{cmd:e(sample)} marks the retained observations after missing values, singleton
chains, and FE/regressor separation are removed. With fweights, {cmd:e(N)} is the
sum of frequency weights rather than the number of marked physical rows.
Weight expressions such as {cmd:[aw=2*w]} are evaluated into temporary doubles.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse ships, clear}{p_end}

{pstd}Basic PPML with fixed effects{p_end}
{phang2}{cmd:. cppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, absorb(ship) vce(robust)}{p_end}

{pstd}With exposure variable{p_end}
{phang2}{cmd:. cppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, absorb(ship) exposure(service) vce(robust)}{p_end}

{pstd}Two-way fixed effects{p_end}
{phang2}{cmd:. cppmlhdfe accident op_75_79, absorb(ship type) vce(cluster ship)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cppmlhdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(rank)}}number of retained coefficients{p_end}
{synopt:{cmd:e(drop_singletons)}}1 if singleton removal was enabled{p_end}
{synopt:{cmd:e(df_a_nested)}}absorbed degrees of freedom nested in clusters{p_end}
{synopt:{cmd:e(df_a_initial)}}absorbed categories before redundancy adjustments{p_end}
{synopt:{cmd:e(df_a_redundant)}}redundant absorbed coefficients{p_end}
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
{synopt:{cmd:e(cmd)}}{cmd:cppmlhdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(absorb)}}absorbed variables{p_end}
{synopt:{cmd:e(d)}}variable containing the saved sum of fixed effects{p_end}
{synopt:{cmd:e(offset)}}offset variable or log-exposure expression{p_end}
{synopt:{cmd:e(exposure)}}exposure variable{p_end}
{synopt:{cmd:e(predict)}}prediction program{p_end}
{synopt:{cmd:e(separation)}}separation methods{p_end}
{synopt:{cmd:e(dofmethod)}}degrees-of-freedom adjustments{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression{p_end}
{synopt:{cmd:e(vce)}}VCE type{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustering){p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(dof_table)}}absorbed categories, redundant coefficients, and nesting by effect{p_end}

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

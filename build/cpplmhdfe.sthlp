{smcl}
{* *! version 1.0.0 11Feb2026}{...}
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
{indepvars}
{ifin}
{cmd:,}
{opt a:bsorb(varlist)}
[{it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(varlist)}}categorical variables representing fixed effects to absorb; required{p_end}
{synopt:{opt exp:osure(varname)}}exposure variable (enters as log offset){p_end}
{synopt:{opt off:set(varname)}}offset variable (enters linearly in the linear predictor){p_end}
{synopt:{opt dof:adjustments(doftype)}}degrees of freedom adjustment method{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt robust} or {opt cluster} {it:clustvar}{p_end}

{syntab:IRLS Convergence}
{synopt:{opt irls_tol:erance(#)}}IRLS convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt irls_max:iter(#)}}maximum IRLS iterations; default is {cmd:1000}{p_end}
{synopt:{opt sep:aration_tolerance(#)}}separation detection tolerance; default is {cmd:1e-8}{p_end}

{syntab:CG Solver}
{synopt:{opt tol:erance(#)}}CG solver convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt iter:ate(#)}}maximum CG iterations; default is {cmd:10000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Reporting}
{synopt:{opt v:erbose}}display detailed timing and convergence information{p_end}
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
are to be absorbed (partialled out). At least one absorb variable is required.

{phang}
{opt exposure(varname)} specifies an exposure variable. This is equivalent to
{opt offset(log(varname))} and is commonly used in rate models.

{phang}
{opt offset(varname)} specifies an offset variable that enters the linear
predictor additively: eta = X*beta + FE + offset.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies the type of standard error. Options are
{opt robust} for Eicker-Huber-White sandwich standard errors and
{opt cluster} {it:clustvar} for cluster-robust standard errors.

{dlgtab:IRLS Convergence}

{phang}
{opt irls_tolerance(#)} specifies the convergence criterion for the IRLS
loop. Convergence is declared when the relative change in deviance falls
below this threshold: |dev - dev_old| / (0.1 + |dev|) < tol. Default is 1e-8.

{phang}
{opt irls_maxiter(#)} specifies the maximum number of IRLS iterations.
Default is 1000.

{phang}
{opt separation_tolerance(#)} specifies the tolerance for detecting
separation. Observations with y=0 and mu < tol are flagged as separated.
Default is 1e-8.


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
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log pseudolikelihood{p_end}
{synopt:{cmd:e(ll_0)}}log pseudolikelihood of null model{p_end}
{synopt:{cmd:e(deviance)}}Poisson deviance{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared (McFadden){p_end}
{synopt:{cmd:e(F)}}Wald F-statistic{p_end}
{synopt:{cmd:e(ic)}}number of IRLS iterations{p_end}
{synopt:{cmd:e(converged)}}1 if IRLS converged, 0 otherwise{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(num_singletons)}}number of singleton observations dropped{p_end}
{synopt:{cmd:e(num_separated)}}number of separated observations dropped{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom absorbed by FEs{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if clustering){p_end}

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

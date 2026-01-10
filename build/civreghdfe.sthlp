{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "civreghdfe##syntax"}{...}
{viewerjumpto "Description" "civreghdfe##description"}{...}
{viewerjumpto "Options" "civreghdfe##options"}{...}
{viewerjumpto "Examples" "civreghdfe##examples"}{...}
{viewerjumpto "Stored results" "civreghdfe##results"}{...}
{title:Title}

{phang}
{bf:civreghdfe} {hline 2} C-accelerated instrumental variables regression with high-dimensional fixed effects


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:civreghdfe}
{it:depvar}
{cmd:(}{it:endogvars} {cmd:=} {it:instruments}{cmd:)}
[{it:exogvars}]
{ifin}
{weight}
{cmd:,}
{opt a:bsorb(varlist)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main (required)}
{synopt:{opt a:bsorb(varlist)}}categorical variables to absorb as fixed effects{p_end}

{syntab:VCE/SE}
{synopt:{opt vce(vcetype)}}variance-covariance estimation: {opt un:adjusted}, {opt r:obust}, {opt cl:uster} {it:clustvar}{p_end}

{syntab:Estimation}
{synopt:{opt tol:erance(#)}}convergence tolerance for CG solver (default: 1e-8){p_end}
{synopt:{opt max:iter(#)}}maximum CG iterations (default: 500){p_end}

{syntab:Reporting}
{synopt:{opt first}}report first-stage regression statistics{p_end}
{synopt:{opt v:erbose}}display progress information{p_end}
{synopt:{opt timeit}}display timing breakdown{p_end}
{synopt:{opt small}}use small-sample adjustments{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{opt aweight}s, {opt fweight}s, and {opt pweight}s are allowed;
see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:civreghdfe} is a high-performance replacement for {cmd:ivreghdfe} that uses
a C plugin for all computations. It performs instrumental variables (IV) regression
with two-stage least squares (2SLS) while absorbing high-dimensional fixed effects.

{pstd}
The command combines the speed of the {cmd:creghdfe} HDFE absorption algorithm
(using conjugate gradient iteration) with 2SLS estimation. All heavy computation
is performed in optimized C code with OpenMP parallelization.

{pstd}
The syntax follows {cmd:ivreghdfe}: endogenous variables and their instruments
are specified in parentheses with an equals sign.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt absorb(varlist)} specifies categorical variables whose fixed effects are
absorbed from all variables (dependent, endogenous, exogenous, and instruments)
before 2SLS estimation. This is required.

{dlgtab:VCE/SE}

{phang}
{opt vce(vcetype)} specifies the variance-covariance estimation method:

{phang2}{opt unadjusted} or {opt ols} - standard VCE assuming homoskedasticity{p_end}
{phang2}{opt robust} - heteroskedasticity-robust standard errors (HC1){p_end}
{phang2}{opt cluster} {it:clustvar} - cluster-robust standard errors{p_end}

{pstd}
The VCE is computed using the proper 2SLS sandwich formula with the original
endogenous regressors (not the first-stage fitted values).

{dlgtab:Estimation}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the conjugate gradient
solver used in HDFE absorption. Default is 1e-8.

{phang}
{opt maxiter(#)} specifies the maximum number of CG iterations. Default is 500.

{dlgtab:Reporting}

{phang}
{opt first} reports first-stage F-statistics for testing instrument strength.

{phang}
{opt verbose} displays detailed progress information during computation.

{phang}
{opt timeit} displays a timing breakdown of computational phases.


{marker examples}{...}
{title:Examples}

{pstd}Basic IV regression with one endogenous variable:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age ttl_exp, absorb(idcode)}{p_end}

{pstd}Multiple endogenous variables:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure hours = union south) age, absorb(idcode year)}{p_end}

{pstd}With robust standard errors:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(robust)}{p_end}

{pstd}With clustered standard errors:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode year) vce(cluster idcode)}{p_end}

{pstd}Display first-stage statistics:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) first}{p_end}

{pstd}With weights:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age [aw=hours], absorb(idcode)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:civreghdfe} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(df_a)}}absorbed degrees of freedom{p_end}
{synopt:{cmd:e(K)}}number of regressors{p_end}
{synopt:{cmd:e(K_endog)}}number of endogenous regressors{p_end}
{synopt:{cmd:e(K_exog)}}number of exogenous regressors{p_end}
{synopt:{cmd:e(K_iv)}}number of instruments (including exogenous){p_end}
{synopt:{cmd:e(G)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if clustered){p_end}
{synopt:{cmd:e(F_first1)}}first-stage F-stat for 1st endogenous var{p_end}
{synopt:{cmd:e(F_first2)}}first-stage F-stat for 2nd endogenous var{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:civreghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(endogvars)}}endogenous variables{p_end}
{synopt:{cmd:e(exogvars)}}exogenous variables{p_end}
{synopt:{cmd:e(instruments)}}excluded instruments{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(vcetype)}}VCE type{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustered){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}


{title:Methods and Formulas}

{pstd}
{cmd:civreghdfe} implements 2SLS estimation with the following steps:

{phang2}1. Absorb fixed effects from all variables (y, endogenous X, exogenous X, instruments Z)
using the conjugate gradient (CG) solver with symmetric Kaczmarz iteration.{p_end}

{phang2}2. Compute the 2SLS estimator on demeaned data:
{it:beta} = (X'P_Z X)^{-1} X'P_Z y
where P_Z = Z(Z'Z)^{-1}Z' is the projection onto the instrument space.{p_end}

{phang2}3. Compute residuals using original X (not projected X) for proper inference.{p_end}

{phang2}4. Compute VCE using the appropriate sandwich estimator for 2SLS.{p_end}

{pstd}
First-stage F-statistics are computed from the partial R-squared of the excluded
instruments in the first-stage regression.


{title:Technical Notes}

{pstd}
The HDFE absorption uses the same algorithm as {cmd:creghdfe}, which iteratively
projects out factor means using a conjugate gradient solver with symmetric
Kaczmarz transformation.

{pstd}
For clustered standard errors, the command uses the standard cluster-robust
sandwich estimator adjusted for 2SLS.


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help ivreghdfe} (if installed), {help creghdfe}, {help ivreg2} (if installed)
{p_end}

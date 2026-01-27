{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "creghdfe##syntax"}{...}
{viewerjumpto "Description" "creghdfe##description"}{...}
{viewerjumpto "Options" "creghdfe##options"}{...}
{viewerjumpto "Examples" "creghdfe##examples"}{...}
{viewerjumpto "Stored results" "creghdfe##results"}{...}
{title:Title}

{phang}
{bf:creghdfe} {hline 2} C-accelerated high-dimensional fixed effects regression

{pstd}
{cmd:creghdfe} is a high-performance replacement for {cmd:reghdfe} version 6.0+
by Sergio Correia ({browse "https://github.com/sergiocorreia/reghdfe":sergiocorreia/reghdfe}).


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:creghdfe}
{depvar}
{indepvars}
{ifin}
{cmd:,}
{opt a:bsorb(varlist)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(varlist[, savefe])}}categorical variables representing fixed effects to absorb; required{p_end}
{synopt:{opt dof:adjustments(doftype)}}degrees of freedom adjustment method{p_end}
{synopt:{opt group:var(newvar)}}save mobility group identifier{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt robust} or {opt cluster} {it:clustvar}{p_end}

{syntab:Convergence}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt max:iter(#)}}maximum iterations; default is {cmd:10000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt nostand:ardize}}do not standardize variables before iteration{p_end}

{syntab:Residuals}
{synopt:{opt resid}}create residual variable named {cmd:_creghdfe_resid}{p_end}
{synopt:{opt resid2(newvar)}}create residual variable with specified name{p_end}
{synopt:{opt resid:uals(newvar)}}alias for {opt resid2()}{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information and timing{p_end}
{synopt:{opt time:it}}display timing breakdown{p_end}
{synoptline}

{pstd}
{it:weight}s are allowed; {opt aweight}s, {opt fweight}s, and {opt pweight}s are supported.


{marker description}{...}
{title:Description}

{pstd}
{cmd:creghdfe} is a high-performance replacement for {help reghdfe:reghdfe}
(if installed) that uses a C plugin with optimized fixed effects absorption.
It estimates linear regression models with multiple high-dimensional fixed effects.

{pstd}
The command uses iterative demeaning to absorb fixed effects, with the
computation performed in C for maximum speed. This is particularly beneficial
for models with many fixed effects or large datasets.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt absorb(varlist[, savefe])} specifies one or more categorical variables whose fixed
effects are to be absorbed. This option is required. The optional {opt savefe}
suboption creates variables {cmd:__hdfe1__}, {cmd:__hdfe2__}, etc. containing
the estimated fixed effect for each observation. {bf:Note:} Variable value storage
is currently a known limitation; the variables are created but may not contain values.

{phang}
{opt dofadjustments(doftype)} controls how degrees of freedom are calculated.
{it:doftype} may be {opt all} (default), {opt none}, {opt firstpair}, or
{opt pairwise}. With {opt none}, no adjustment for connected groups is made.
With {opt firstpair}, only the first two FE variables are checked for
connectivity. With {opt pairwise} or {opt all}, all FE pairs are checked.

{phang}
{opt groupvar(newvar)} creates a new variable containing the mobility group
identifier for each observation. Mobility groups are connected components
in the network of fixed effects. {bf:Note:} Variable value storage is currently
a known limitation; the variable is created but may not contain values.

{dlgtab:SE/Robust}

{phang}
{opt vce(robust)} computes heteroskedasticity-robust standard errors.

{phang}
{opt vce(cluster clustvar)} computes cluster-robust standard errors, clustering
on the specified variable.

{dlgtab:Convergence}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the iterative
demeaning algorithm. The default is {cmd:1e-8}.

{phang}
{opt maxiter(#)} specifies the maximum number of iterations for the
demeaning algorithm. The default is {cmd:10000}.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations including fixed effects absorption and OLS computation. By default,
{cmd:creghdfe} uses all available CPU cores as reported by OpenMP. Use this
option to limit parallelism, for example when running multiple jobs simultaneously.

{phang}
{opt nostandardize} specifies that variables should not be standardized
before the iterative algorithm. This may affect convergence speed.

{dlgtab:Residuals}

{phang}
{opt resid} stores residuals in a new variable named {cmd:_creghdfe_resid}.

{phang}
{opt resid2(newvar)} stores residuals in a new variable with the specified name.

{phang}
{opt residuals(newvar)} is an alias for {opt resid2(newvar)}, provided for
compatibility with {help reghdfe:reghdfe}.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information during estimation,
including timing for each stage of the computation.

{phang}
{opt timeit} displays timing breakdown for the estimation.


{marker examples}{...}
{title:Examples}

{pstd}Setup:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Basic regression with a single fixed effect:{p_end}
{phang2}{cmd:. creghdfe price mpg weight, absorb(foreign)}{p_end}

{pstd}With robust standard errors:{p_end}
{phang2}{cmd:. creghdfe price mpg weight, absorb(foreign) vce(robust)}{p_end}

{pstd}Using panel data with two-way fixed effects:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year)}{p_end}

{pstd}With clustered standard errors:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)}{p_end}

{pstd}Save residuals:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp, absorb(idcode) resid2(myresid)}{p_end}

{pstd}With verbose output showing timing:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:creghdfe} stores results in {cmd:e()} similar to {help reghdfe:reghdfe}.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command provides similar functionality to {cmd:reghdfe} by Sergio Correia.
The original {cmd:reghdfe} package is available at
{browse "https://github.com/sergiocorreia/reghdfe"}. We thank Sergio Correia
for his groundbreaking work on high-dimensional fixed effects estimation in Stata.

{pstd}
See also: Correia, S. 2016. "Linear Models with High-Dimensional Fixed Effects:
An Efficient and Feasible Estimator." Working paper.
{browse "http://scorreia.com/research/hdfe.pdf"}


{title:Also see}

{psee}
Manual: {bf:[R] areg}

{psee}
Online: {browse "https://github.com/sergiocorreia/reghdfe":reghdfe} (if installed),
{help areg}, {help xtreg}, {help ctools}
{p_end}

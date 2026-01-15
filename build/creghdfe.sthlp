{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "creghdfe##syntax"}{...}
{viewerjumpto "Description" "creghdfe##description"}{...}
{viewerjumpto "Options" "creghdfe##options"}{...}
{viewerjumpto "Examples" "creghdfe##examples"}{...}
{viewerjumpto "Stored results" "creghdfe##results"}{...}
{title:Title}

{phang}
{bf:creghdfe} {hline 2} C-accelerated high-dimensional fixed effects regression


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
{synopt:{opt a:bsorb(varlist)}}categorical variables representing fixed effects to absorb; required{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt robust} or {opt cluster} {it:clustvar}{p_end}

{syntab:Convergence}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt max:iter(#)}}maximum iterations; default is {cmd:10000}{p_end}
{synopt:{opt nostand:ardize}}do not standardize variables before iteration{p_end}

{syntab:Residuals}
{synopt:{opt resid}}create residual variable named {cmd:_creghdfe_resid}{p_end}
{synopt:{opt resid2(newvar)}}create residual variable with specified name{p_end}

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
{opt absorb(varlist)} specifies one or more categorical variables whose fixed
effects are to be absorbed. This option is required.

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
{opt nostandardize} specifies that variables should not be standardized
before the iterative algorithm. This may affect convergence speed.

{dlgtab:Residuals}

{phang}
{opt resid} stores residuals in a new variable named {cmd:_creghdfe_resid}.

{phang}
{opt resid2(newvar)} stores residuals in a new variable with the specified name.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information during estimation,
including timing for each stage of the computation.

{phang}
{opt timeit} displays timing breakdown for the estimation.


{marker examples}{...}
{title:Examples}

{pstd}Basic regression with firm and year fixed effects:{p_end}
{phang2}{cmd:. creghdfe y x1 x2, absorb(firm year)}{p_end}

{pstd}With clustered standard errors:{p_end}
{phang2}{cmd:. creghdfe y x1 x2, absorb(firm year) vce(cluster firm)}{p_end}

{pstd}With verbose output:{p_end}
{phang2}{cmd:. creghdfe y x1 x2, absorb(firm year) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:creghdfe} stores results in {cmd:e()} similar to {help reghdfe:reghdfe}.


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help reghdfe} (if installed), {help ctools}, {help areg}
{p_end}

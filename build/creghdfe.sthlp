{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "creghdfe##syntax"}{...}
{viewerjumpto "Description" "creghdfe##description"}{...}
{viewerjumpto "Options" "creghdfe##options"}{...}
{viewerjumpto "Prediction" "creghdfe##prediction"}{...}
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
{weight}
[{cmd:,} {opt a:bsorb(varlist)} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(varlist[, savefe])}}categorical variables representing fixed effects to absorb; optional{p_end}
{synopt:{opt dof:adjustments(doftype)}}degrees of freedom adjustment method{p_end}
{synopt:{opt group:var(newvar)}}save mobility group identifier{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt unadjusted} (default), {opt robust}, or {opt cluster} {it:clustvar}{p_end}

{syntab:Convergence}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt iterate(#)}}maximum iterations; default is {cmd:10000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt nostand:ardize}}do not standardize variables before iteration{p_end}

{syntab:Residuals}
{synopt:{opt resid}}create residual variable named {cmd:_reghdfe_resid}{p_end}
{synopt:{opt resid2(newvar)}}create residual variable with specified name{p_end}
{synopt:{opt resid:uals(newvar)}}alias for {opt resid2()}{p_end}

{syntab:Advanced}
{synopt:{opt quad}}use quad-precision accumulation in CG solver{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information and timing{p_end}
{synoptline}

{pstd}
{it:weight}s are allowed; {opt aweight}s, {opt fweight}s, and {opt pweight}s are supported.
{opt pweight}s imply {opt vce(robust)} unless clustering is specified.


{marker description}{...}
{title:Description}

{pstd}
{cmd:creghdfe} estimates linear regressions with multiple categorical fixed
effects. The fixed effects are removed by iterative projection in C. Without
{opt absorb()}, the command estimates a linear regression with a constant.
Factor-variable and time-series notation are allowed for the regressors.

{pstd}
Observations with missing model variables are excluded. Singleton fixed-effect
groups are dropped repeatedly until none remain. Collinear regressors are
omitted. {cmd:e(sample)} marks the observations retained for estimation.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt absorb(varlist[, savefe])} specifies one or more categorical variables whose fixed
effects are to be absorbed. Each term
is a numeric or string variable, an {cmd:i.}{it:var}, or an interaction
{it:var1}{cmd:#}{it:var2} absorbed as one effect; heterogeneous slopes
({cmd:c.}{it:var}), {cmd:##} and {it:newvar}{cmd:=} terms are not supported. The optional {opt savefe}
suboption creates variables {cmd:__hdfe1__}, {cmd:__hdfe2__}, etc. containing
the estimated fixed effect for each retained observation. The output names
must be new. Individual fixed-effect levels depend on normalization.

{phang}
{opt dofadjustments(list)} selects how redundant fixed-effect coefficients are
counted, as in {cmd:reghdfe}. {it:list} may be {opt all} (default), {opt none},
or any of {opt first:pair} or {opt pair:wise}, {opt cl:usters} and
{opt cont:inuous}. {opt pairwise} counts each effect's redundant coefficients
from the mobility groups it forms with every earlier effect; {opt firstpair}
examines the first two effects only. {opt clusters} treats effects nested
within the cluster variable as costing no degrees of freedom and sets them
aside before mobility groups are counted. {opt all} is {opt pairwise clusters}
{opt continuous}. {opt continuous} concerns slope coefficients, which
{cmd:creghdfe} does not absorb, and has no effect.

{phang}
{opt groupvar(newvar)} creates a new variable containing the mobility group
identifier for each observation. Mobility groups are connected components in the network formed by the first
two fixed effects. With one effect, all retained observations receive 1.

{dlgtab:SE/Robust}

{phang}
{opt vce(robust)} computes heteroskedasticity-robust standard errors.

{phang}
{opt vce(cluster clustvar)} computes cluster-robust standard errors, clustering
on the specified variable, which may be a string variable, an {cmd:i.}{it:var}, or
an interaction {it:var1}{cmd:#}{it:var2}. Only one cluster dimension is allowed;
multi-way clustering is not supported.

{phang}
With {cmd:pweight}s the standard errors are robust unless {opt vce(cluster)} is
specified; {opt vce(unadjusted)} and {opt vce(ols)} are an error, as in {cmd:reghdfe}.

{dlgtab:Convergence}

{phang}
{opt tolerance(#)} specifies a positive convergence tolerance for the iterative
demeaning algorithm. The default is {cmd:1e-8}.

{phang}
{opt iterate(#)} specifies a positive maximum number of iterations for the
demeaning algorithm. The default is {cmd:10000}.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt nostandardize} specifies that variables should not be standardized
before the iterative algorithm. This may affect convergence speed.

{dlgtab:Advanced}

{phang}
{opt quad} uses quad-precision (double-double) accumulation in the conjugate
gradient solver for improved numerical accuracy. This is slower but may help
with ill-conditioned problems.

{dlgtab:Residuals}

{phang}
{opt resid} stores residuals in a new variable named {cmd:_reghdfe_resid}.

{phang}
{opt resid2(newvar)} stores residuals in a new variable with the specified name.

{phang}
{opt residuals(newvar)} is an alias for {opt resid2(newvar)}, provided for
compatibility with {help reghdfe:reghdfe}.

{dlgtab:Data recovery}

{pstd}
By default, an error or interruption removes new output variables and clears
{cmd:e()}. An exhausted projection solver returns error 430.

{phang}
{opt snapshot} uses {help preserve} to copy the whole dataset for recovery.
Stata may write the copy to disk when it exceeds {cmd:c(max_preservemem)}.

{phang}
{opt fast} skips data recovery. An error or interruption can leave partially
filled output variables. It cannot be combined with {opt snapshot}.


{dlgtab:Reporting}

{phang}
{opt verbose} displays progress and elapsed times.


{marker prediction}{...}
{title:Prediction}

{pstd}
{cmd:predict newvar} defaults to {opt xb}, the linear prediction excluding
absorbed effects. {opt stdp} gives its standard error. Estimate with {opt resid},
{opt resid2()}, or {opt residuals()} to use the following prediction options:

{p 8 12 2}{opt residuals}: stored regression residuals.{p_end}
{p 8 12 2}{opt d}: sum of absorbed effects, y - xb - residual.{p_end}
{p 8 12 2}{opt xbd}: fitted values including absorbed effects, y - residual.{p_end}
{p 8 12 2}{opt dresiduals}: absorbed effects plus residual, y - xb.{p_end}

{pstd}
Residual and fixed-effect output variables are missing outside the estimation
sample. Their names, including those from {opt groupvar()} and {opt savefe},
must be distinct new variables.

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
{cmd:creghdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_within)}}within R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(tss_within)}}within total sum of squares{p_end}
{synopt:{cmd:e(df_a)}}absorbed degrees of freedom (adjusted){p_end}
{synopt:{cmd:e(df_a_initial)}}levels of the absorbed effects, before redundant coefficients are removed{p_end}
{synopt:{cmd:e(df_a_redundant)}}redundant coefficients of the absorbed effects, nested ones included{p_end}
{synopt:{cmd:e(df_a_nested)}}degrees of freedom from nested FEs{p_end}
{synopt:{cmd:e(F)}}F-statistic{p_end}
{synopt:{cmd:e(rank)}}number of retained slope coefficients{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(N_hdfe_extended)}}number of absorbed effects, as in {cmd:reghdfe} (the mobility groups are the second row's redundant count in {cmd:e(dof_table)}){p_end}
{synopt:{cmd:e(num_singletons)}}number of singleton observations dropped{p_end}
{synopt:{cmd:e(drop_singletons)}}1: singletons are always dropped{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(N_clust1)}}number of clusters (if {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster dimensions (if {cmd:vce(cluster)}){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:creghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effect variables{p_end}
{synopt:{cmd:e(absvars)}}absorbed effects as specified ({cmd:_cons} without {opt absorb()}){p_end}
{synopt:{cmd:e(extended_absvars)}}absorbed effects without {cmd:i.} prefixes{p_end}
{synopt:{cmd:e(title)}}{cmd:HDFE Linear regression}{p_end}
{synopt:{cmd:e(title2)}}number of absorbed effects{p_end}
{synopt:{cmd:e(dofmethod)}}degrees-of-freedom adjustments applied{p_end}
{synopt:{cmd:e(vce)}}VCE type: {cmd:unadjusted}, {cmd:robust}, or {cmd:cluster}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust} (if {cmd:vce(robust)} or {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable (if {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(clustvar1)}}name of cluster variable (if {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(resid)}}name of residual variable (if {cmd:resid} or {cmd:resid2()} specified){p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(dof_table)}}per effect: categories, redundant, coefficients, inexact (missing: not computed) and nested{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command provides similar functionality to {cmd:reghdfe} by Sergio Correia.
The original {cmd:reghdfe} package is available at
{browse "https://github.com/sergiocorreia/reghdfe"}.

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

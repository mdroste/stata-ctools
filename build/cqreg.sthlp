{smcl}
{* *! version 1.2.4 20261103}{...}
{vieweralsosee "[R] qreg" "help qreg"}{...}
{vieweralsosee "ctools" "help ctools"}{...}
{viewerjumpto "Syntax" "cqreg##syntax"}{...}
{viewerjumpto "Description" "cqreg##description"}{...}
{viewerjumpto "Options" "cqreg##options"}{...}
{viewerjumpto "Prediction" "cqreg##prediction"}{...}
{viewerjumpto "Examples" "cqreg##examples"}{...}
{viewerjumpto "Stored results" "cqreg##results"}{...}
{viewerjumpto "Methods" "cqreg##methods"}{...}
{title:Title}

{phang}
{bf:cqreg} {hline 2} C-accelerated quantile regression


{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:cqreg}
{depvar}
[{indepvars}]
{ifin}
[{cmd:,} {it:options}]


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt q:uantile(#)}}estimate {it:#} quantile; default is {cmd:quantile(0.5)}{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimation method; {it:vcetype} may be
    {opt iid}, {opt robust}, or {opt cl:uster} {it:clustvar}{p_end}
{synopt:{opt den:method(method)}}density estimation: {opt fitted}, {opt residual}, or {opt kernel}[{cmd:(}{it:kernel}{cmd:)}]{p_end}
{synopt:{opt bwmethod(method)}}bandwidth selection: {opt hsheather}, {opt bofinger}, or {opt chamberlain}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{it:display_options}}{help estimation options##display_options:display options} of {cmd:qreg}{p_end}
{synopt:{opt nolog}}accepted for compatibility with {cmd:qreg}; {cmd:cqreg} prints no iteration log{p_end}
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Optimization}
{synopt:{opt tol:erance(#)}}convergence tolerance in normalized response units; default is {cmd:tolerance(1e-12)}{p_end}
{synopt:{opt max:iter(#)}}maximum iterations; default is {cmd:maxiter(200)}{p_end}
{synopt:{opt nopreprocess(#)}}accepted for backward compatibility; has no effect{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cqreg} estimates a conditional quantile of the dependent variable by
minimizing the quantile loss. The default is median regression. It uses a
primal-dual interior-point solver in C.

{pstd}
Factor-variable notation is allowed. Include explicit indicators for fixed
effects, for example {cmd:cqreg y x i.group}. Without {it:indepvars},
{cmd:cqreg} fits the constant-only model. Weights and {opt absorb()} are
not accepted. Typing {cmd:cqreg} without arguments, or
{cmd:estimates replay}, redisplays the last estimates; {opt level()} and the
display options may be given then.

{pstd}
The quantile regression objective can have several minimizers. With indicator
regressors this is common: when a group's size times {it:q} is an integer,
that group's coefficient can move within an interval without changing the
objective. When the minimum is not unique, {cmd:cqreg} and {cmd:qreg} may
report different coefficients with the same quantile loss, stored in
{cmd:e(sum_adev)}. {cmd:cqreg} checks for alternative solutions and reports
those it finds; see {help cqreg##methods:Methods and formulas} for the limits
of this check.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt quantile(#)} specifies the quantile to be estimated. The quantile must be
a number strictly between 0 and 1. The default is {cmd:quantile(0.5)}, which
corresponds to the median.


{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies how standard errors are computed.

{phang2}
{opt vce(iid)} (the default) computes standard errors assuming independent and
identically distributed errors. The sparsity function is estimated by the
method in {opt denmethod()}.

{phang2}
{opt vce(robust)} computes heteroskedasticity-robust standard errors using the
Powell sandwich estimator.

{phang2}
{cmd:vce(cluster} {it:clustvar}{cmd:)} computes cluster-robust standard errors
clustered on {it:clustvar}. With a single cluster the covariance is undefined:
as {cmd:qreg}, {cmd:cqreg} posts a zero {cmd:e(V)} with {cmd:e(rank)} 0 and
displays no standard errors.

{phang}
{opt denmethod(method)} specifies the density estimation method for computing
standard errors; it corresponds to the suboptions of {help qreg}'s
{cmd:vce(}{it:vcetype}{cmd:,} {it:method}{cmd:)}.

{phang2}
{opt denmethod(fitted)} (the default) uses the fitted (Siddiqui) method based on
the difference quotient of the quantile function.

{phang2}
{opt denmethod(residual)} uses the difference quotient of the percentiles of the
residuals, excluding the observations of the linear programming basis, as
{cmd:qreg}'s {cmd:vce(iid, residual)}. It may be combined only with
{cmd:vce(iid)}.

{phang2}
{cmd:denmethod(kernel}[{cmd:(}{it:kernel}{cmd:)}]{cmd:)} uses the kernel density
estimate of the residuals, as {cmd:qreg}'s {cmd:vce(}{it:vcetype}{cmd:, kernel(}{it:kernel}{cmd:))}.
{it:kernel} is {opt ep:anechnikov} (the default), {opt epan2}, {opt bi:weight},
{opt cos:ine}, {opt gau:ssian}, {opt par:zen}, {opt rec:tangle}, or
{opt tri:angle}.

{phang}
{opt bwmethod(method)} specifies the bandwidth selection method for the
density estimation used in computing standard errors.

{phang2}
{opt bwmethod(hsheather)} (the default) uses the Hall-Sheather (1988) bandwidth.

{phang2}
{opt bwmethod(bofinger)} uses the Bofinger (1975) bandwidth.

{phang2}
{opt bwmethod(chamberlain)} uses the Chamberlain (1994) bandwidth.

{dlgtab:Reporting}

{phang}
{opt level(#)} sets the confidence level of the reported intervals. As in
{cmd:qreg}, it also enters the Hall-Sheather and Chamberlain bandwidths, and
so the standard errors, through alpha = 1 - {it:#}/100.

{phang}
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{dlgtab:Optimization}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the IPM solver. The
default is {cmd:tolerance(1e-12)}. The response is centered at its median and
scaled by a power of two before fitting, so its median absolute deviation is
between 1 and 2 (the mean absolute deviation when most values are tied). Data
whose spread is already in that range are only centered. Centering reduces
rounding error when the dependent variable has a large level relative to its
variation. Convergence is assessed in these normalized units. Coefficients,
covariance, sparsity and deviation statistics are returned in the original
response units. If covariance values are too small or too large to calculate
or store in Stata, {cmd:cqreg} returns error 498 without storing new
estimates. Rescale the dependent variable and estimate the model again.

{phang}
{opt maxiter(#)} sets the maximum number of iterations of the interior-point
algorithm. The default is 200. If the algorithm fails to converge or cannot
calculate finite estimates, {cmd:cqreg} returns error 430 without storing new
estimates.

{phang}
{opt nopreprocess(#)} accepts nonnegative integers for compatibility and has
no effect on the solver.

{marker prediction}{...}
{title:Prediction}

{pstd}
{cmd:predict newvar} defaults to {opt xb}, the fitted conditional quantile.
{opt residuals} returns y - xb, and {opt stdp} returns the standard error of
the linear prediction. {cmd:if} and {cmd:in} may restrict the predictions.


{marker examples}{...}
{title:Examples}

{pstd}Setup:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Basic median regression:{p_end}
{phang2}{cmd:. cqreg price mpg weight}{p_end}

{pstd}25th percentile regression:{p_end}
{phang2}{cmd:. cqreg price mpg weight, quantile(0.25)}{p_end}

{pstd}75th percentile with robust standard errors:{p_end}
{phang2}{cmd:. cqreg price mpg weight, quantile(0.75) vce(robust)}{p_end}

{pstd}Robust standard errors from a biweight kernel density estimate:{p_end}
{phang2}{cmd:. cqreg price mpg weight, vce(robust) denmethod(kernel(biweight))}{p_end}

{pstd}Quantile regression with fixed effects:{p_end}
{phang2}{cmd:. webuse grunfeld, clear}{p_end}
{phang2}{cmd:. cqreg invest mvalue kstock i.company}{p_end}

{pstd}Two-way fixed effects with clustering:{p_end}
{phang2}{cmd:. cqreg invest mvalue kstock i.company i.year, vce(cluster company)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cqreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(q)}}quantile{p_end}
{synopt:{cmd:e(sum_adev)}}minimized quantile loss, sum of q*u + (1-q)*v{p_end}
{synopt:{cmd:e(sum_rdev)}}quantile loss from the intercept-only model{p_end}
{synopt:{cmd:e(r2_p)}}1 - e(sum_adev)/e(sum_rdev){p_end}
{synopt:{cmd:e(q_v)}}unconditional sample quantile of the dependent variable{p_end}
{synopt:{cmd:e(f_r)}}estimated density, 1/e(sparsity){p_end}
{synopt:{cmd:e(sum_w)}}number of observations (weights are not accepted){p_end}
{synopt:{cmd:e(sparsity)}}estimated sparsity (1/f(0)); missing with
{cmd:vce(robust)} or {cmd:vce(cluster)}, as in {cmd:qreg}; zero when the
dependent variable has a mass point at the fitted quantile, as in {cmd:qreg},
which then gives a zero covariance{p_end}
{synopt:{cmd:e(bwidth)}}bandwidth used for sparsity estimation{p_end}
{synopt:{cmd:e(kbwidth)}}kernel bandwidth (kernel method){p_end}
{synopt:{cmd:e(iterations)}}number of IPM iterations{p_end}
{synopt:{cmd:e(convcode)}}0 for a successfully converged fit{p_end}
{synopt:{cmd:e(nonunique)}}1 if the check finds other coefficients with the same minimum;
0 if it does not; missing if the check could not be performed. A value of 0
does not establish that the solution is unique{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom, N - {cmd:e(rank)}{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}; 0 for a zero covariance (a single
cluster, or zero sparsity), as in {cmd:qreg}{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if {opt vce(cluster)}){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cqreg}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(vce)}}variance estimation method{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust} with {cmd:vce(robust)} or {cmd:vce(cluster)}{p_end}
{synopt:{cmd:e(denmethod)}}density estimation method{p_end}
{synopt:{cmd:e(kernel)}}kernel function (kernel method){p_end}
{synopt:{cmd:e(bwmethod)}}bandwidth selection method{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if {opt vce(cluster)}){p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{cmd:cqreg} minimizes the quantile regression objective function:

{phang2}
minimize {it:sum_i} [ {it:q} * {it:u_i} + (1-{it:q}) * {it:v_i} ]

{phang2}
subject to: {it:y_i} - {it:x_i}'{it:beta} = {it:u_i} - {it:v_i}, {it:u_i} >= 0, {it:v_i} >= 0

{pstd}
where {it:q} is the quantile, {it:u_i} represents positive deviations, and
{it:v_i} represents negative deviations.

{pstd}
The solver uses primal-dual iterations with Mehrotra predictor-corrector steps,
then moves to a vertex of the linear program: the fit through K observations
(the basis) with zero residuals.

{pstd}
The optimum need not be unique, which is common with indicator regressors (a
group whose size times q is an integer). At the vertex, each edge of the
linear program moves one basis residual off zero; {cmd:cqreg} computes the
rate at which the objective changes along each edge (the reduced costs) and,
when a rate is zero to rounding, sets {cmd:e(nonunique)} to 1 and notes that
other coefficients attain the same minimum. {cmd:qreg} and {cmd:cqreg} can
then report different, equally optimal coefficients. The check can miss
alternative solutions when additional observations have zero residuals. A
finding of nonuniqueness identifies an alternative solution; failure to find
one does not establish uniqueness.

{pstd}
Standard errors follow {cmd:qreg}. The IID variance formula is:

{phang2}
V = sparsity^2 * q*(1-q) * (X'X)^{-1}

{pstd}
where the sparsity 1/f(0) is estimated with bandwidth h ({opt bwmethod()}) by
the method in {opt denmethod()}: {opt fitted} evaluates quantile regressions at
q-h and q+h at the mean of X; {opt residual} takes the difference quotient
(F^{-1}(q+h) - F^{-1}(q-h))/(2h) of the residual percentiles after dropping the
K observations of the linear programming basis (h uses n-K observations);
{opt kernel} uses kband/mean(K(r_i/kband)), where
kband = min(sd, IQR/1.34) * (invnormal(q+h) - invnormal(q-h)) of the residuals r_i.
The robust VCE is q*(1-q) * H * (X'X) * H, with H the inverse of
sum_i f_i x_i x_i' ({opt fitted}) or of sum_i K(r_i/kband) x_i x_i' / kband
({opt kernel}).

{pstd}
Fixed effects must enter as explicit indicators in the quantile objective;
least-squares partialling is not used.


{marker references}{...}
{title:References}

{phang}
Koenker, R. 2005. {it:Quantile Regression}. Cambridge University Press.

{phang}
Chernozhukov, V., I. Fernández-Val, and B. Melly. 2020. Fast algorithms for
the quantile regression process. {it:Empirical Economics} 62: 7-33.

{phang}
Portnoy, S., and R. Koenker. 1997. The Gaussian hare and the Laplacian tortoise:
computability of squared-error versus absolute-error estimators.
{it:Statistical Science} 12: 279-300.

{phang}
Hall, P., and S. J. Sheather. 1988. On the distribution of a studentized quantile.
{it:Journal of the Royal Statistical Society, Series B} 50: 381-391.


{marker author}{...}
{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{marker seealso}{...}
{title:Also see}

{psee}
Manual: {bf:[R] qreg}

{psee}
Online: {help qreg}, {help bsqreg}, {help sqreg}, {help iqreg}, {help ctools}
{p_end}

{smcl}
{* *! version 1.1.0 20261016}{...}
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
{indepvars}
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
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Optimization}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:tolerance(1e-12)}{p_end}
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
effects, for example {cmd:cqreg y x i.group}. Weights and {opt absorb()} are
not accepted.


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
clustered on {it:clustvar}.

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
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{dlgtab:Optimization}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the IPM solver.
The default is {cmd:tolerance(1e-12)}.

{phang}
{opt maxiter(#)} specifies the maximum number of IPM iterations. The default
is {cmd:maxiter(200)}. If the solver exhausts this limit, encounters a numerical
breakdown, or produces nonfinite results, {cmd:cqreg} returns error 430 and does
not post new estimates.

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
{synopt:{cmd:e(sparsity)}}estimated sparsity (1/f(0)); missing for the kernel
method with {cmd:vce(robust)} or {cmd:vce(cluster)}, as in {cmd:qreg}{p_end}
{synopt:{cmd:e(bwidth)}}bandwidth used for sparsity estimation{p_end}
{synopt:{cmd:e(kbwidth)}}kernel bandwidth (kernel method){p_end}
{synopt:{cmd:e(iterations)}}number of IPM iterations{p_end}
{synopt:{cmd:e(convcode)}}0 for a successfully converged fit{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if {opt vce(cluster)}){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cqreg}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(vce)}}variance estimation method{p_end}
{synopt:{cmd:e(denmethod)}}density estimation method{p_end}
{synopt:{cmd:e(kernel)}}kernel function (kernel method){p_end}
{synopt:{cmd:e(bwmethod)}}bandwidth selection method{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if {opt vce(cluster)}){p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}

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
The solver uses primal-dual iterations with Mehrotra predictor-corrector steps.

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

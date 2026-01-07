{smcl}
{* *! version 1.0.0 07Jan2026}{...}
{vieweralsosee "[R] qreg" "help qreg"}{...}
{vieweralsosee "ctools" "help ctools"}{...}
{viewerjumpto "Syntax" "cqreg##syntax"}{...}
{viewerjumpto "Description" "cqreg##description"}{...}
{viewerjumpto "Options" "cqreg##options"}{...}
{viewerjumpto "Examples" "cqreg##examples"}{...}
{viewerjumpto "Stored results" "cqreg##results"}{...}
{viewerjumpto "Methods" "cqreg##methods"}{...}
{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{bf:cqreg} {hline 2}}C-accelerated quantile regression{p_end}
{p2colreset}{...}


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
{synopt:{opt abs:orb(varlist)}}categorical variables to absorb (HDFE){p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimation method; {it:vcetype} may be
    {opt iid}, {opt robust}, or {opt cl:uster} {it:clustvar}{p_end}
{synopt:{opt bwmethod(method)}}bandwidth selection: {opt hsheather}, {opt bofinger}, or {opt chamberlain}{p_end}

{syntab:Reporting}
{synopt:{opt v:erbose}}display progress information{p_end}
{synopt:{opt timeit}}display execution time{p_end}

{syntab:Optimization}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:tolerance(1e-8)}{p_end}
{synopt:{opt max:iter(#)}}maximum iterations; default is {cmd:maxiter(50)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cqreg} fits quantile regression models using a high-performance C implementation
with an Interior Point Method (IPM) solver. It is a drop-in replacement for Stata's
{helpb qreg} command with significant performance improvements, especially for
large datasets.

{pstd}
Key features:

{phang2}1. {bf:Interior Point Method solver} - Primal-dual IPM with Mehrotra predictor-corrector
for fast convergence, especially on large datasets (N > 100,000).{p_end}

{phang2}2. {bf:HDFE support} - Unlike native {cmd:qreg}, {cmd:cqreg} can absorb high-dimensional
fixed effects via the {opt absorb()} option using a conjugate gradient solver.{p_end}

{phang2}3. {bf:Parallel computation} - Uses OpenMP for parallel linear algebra operations
when available.{p_end}

{phang2}4. {bf:Multiple VCE options} - Supports IID (default), robust, and clustered
standard errors.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt quantile(#)} specifies the quantile to be estimated. The quantile must be
a number strictly between 0 and 1. The default is {cmd:quantile(0.5)}, which
corresponds to the median.

{phang}
{opt absorb(varlist)} specifies categorical variables whose fixed effects are
to be absorbed (partialled out) from the regression. This enables estimation
of quantile regression with high-dimensional fixed effects without creating
dummy variables. The variables should be integer-valued group identifiers.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies how standard errors are computed.

{phang2}
{opt vce(iid)} (the default) computes standard errors assuming independent and
identically distributed errors. Uses kernel density estimation for the sparsity
function.

{phang2}
{opt vce(robust)} computes heteroskedasticity-robust standard errors using the
Powell sandwich estimator.

{phang2}
{opt vce(cluster} {it:clustvar}{opt )} computes cluster-robust standard errors
clustered on {it:clustvar}.

{phang}
{opt bwmethod(method)} specifies the bandwidth selection method for kernel
density estimation used in computing standard errors.

{phang2}
{opt bwmethod(hsheather)} (the default) uses the Hall-Sheather (1988) bandwidth.

{phang2}
{opt bwmethod(bofinger)} uses the Bofinger (1975) bandwidth.

{phang2}
{opt bwmethod(chamberlain)} uses the Chamberlain (1994) bandwidth.

{dlgtab:Reporting}

{phang}
{opt verbose} displays progress information during estimation, including
timing breakdowns.

{phang}
{opt timeit} displays the total execution time.

{dlgtab:Optimization}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the IPM solver.
The default is {cmd:tolerance(1e-8)}.

{phang}
{opt maxiter(#)} specifies the maximum number of IPM iterations. The default
is {cmd:maxiter(50)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Basic median regression{p_end}
{phang2}{cmd:. cqreg price mpg weight}{p_end}

{pstd}25th percentile regression{p_end}
{phang2}{cmd:. cqreg price mpg weight, quantile(0.25)}{p_end}

{pstd}75th percentile with robust standard errors{p_end}
{phang2}{cmd:. cqreg price mpg weight, quantile(0.75) vce(robust)}{p_end}

{pstd}Quantile regression with fixed effects{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. cqreg ln_wage age ttl_exp tenure, absorb(idcode)}{p_end}

{pstd}Two-way fixed effects with clustering{p_end}
{phang2}{cmd:. cqreg ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cqreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(q)}}quantile{p_end}
{synopt:{cmd:e(sum_adev)}}sum of absolute deviations{p_end}
{synopt:{cmd:e(sparsity)}}estimated sparsity (1/f(0)){p_end}
{synopt:{cmd:e(bwidth)}}bandwidth used for sparsity estimation{p_end}
{synopt:{cmd:e(iterations)}}number of IPM iterations{p_end}
{synopt:{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom absorbed by FEs (if {opt absorb()}){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if {opt vce(cluster)}){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cqreg}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(vce)}}variance estimation method{p_end}
{synopt:{cmd:e(bwmethod)}}bandwidth selection method{p_end}
{synopt:{cmd:e(absorb)}}absorbed FE variables (if specified){p_end}
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
The Interior Point Method solves this linear program using a primal-dual
approach with Mehrotra predictor-corrector steps. This provides O(sqrt(N))
convergence, which is significantly faster than the simplex method used by
native {cmd:qreg} for large datasets.

{pstd}
Standard errors are computed using kernel density estimation for the sparsity
function f(0) = density at the conditional quantile. The IID variance formula is:

{phang2}
V = (1/n) * sparsity^2 * q*(1-q) * (X'X)^{-1}

{pstd}
For HDFE models, the fixed effects are partialled out using a conjugate gradient
solver with symmetric Kaczmarz transformations before applying the IPM solver.


{marker references}{...}
{title:References}

{phang}
Koenker, R. 2005. {it:Quantile Regression}. Cambridge University Press.

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
Part of the ctools suite.


{marker seealso}{...}
{title:Also see}

{psee}
Manual: {bf:[R] qreg}

{psee}
{space 2}Help: {helpb qreg}, {helpb bsqreg}, {helpb sqreg}, {helpb iqreg}
{p_end}

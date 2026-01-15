{smcl}
{* *! version 1.1.0}{...}
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

{syntab:Estimators}
{synopt:{opt liml}}use limited-information maximum likelihood (LIML){p_end}
{synopt:{opt fuller(#)}}use Fuller's modified LIML with parameter {it:#}{p_end}
{synopt:{opt kclass(#)}}use k-class estimator with {it:k} = {it:#}{p_end}
{synopt:{opt gmm2s}}use two-step efficient GMM{p_end}
{synopt:{opt cue}}use continuously-updated GMM (CUE){p_end}

{syntab:VCE/SE}
{synopt:{opt vce(vcetype)}}variance-covariance estimation: {opt un:adjusted}, {opt r:obust}, {opt cl:uster} {it:clustvar}{p_end}
{synopt:{opt small}}use small-sample adjustments (df corrections){p_end}

{syntab:HAC Standard Errors}
{synopt:{opt bw(#)}}bandwidth for kernel-based HAC estimation{p_end}
{synopt:{opt kernel(string)}}kernel type for HAC: {opt bartlett}, {opt parzen}, {opt quadraticspectral}, {opt truncated}, {opt tukey}{p_end}
{synopt:{opt dkraay(#)}}Driscoll-Kraay SEs with {it:#} lags (for panel data){p_end}

{syntab:Estimation Settings}
{synopt:{opt tol:erance(#)}}convergence tolerance for CG solver (default: 1e-8){p_end}
{synopt:{opt max:iter(#)}}maximum CG iterations (default: 500){p_end}
{synopt:{opt noc:onstant}}suppress constant term (absorbed with FE){p_end}

{syntab:Reporting}
{synopt:{opt first}}report first-stage regression statistics{p_end}
{synopt:{opt rf}}report reduced-form estimates{p_end}
{synopt:{opt coviv}}display covariance matrix of IV estimators{p_end}
{synopt:{opt v:erbose}}display progress information{p_end}
{synopt:{opt timeit}}display timing breakdown{p_end}

{syntab:Save Results}
{synopt:{opt res:iduals(newvar)}}save residuals to new variable{p_end}
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
using 2SLS, LIML, k-class, GMM, or CUE estimators while absorbing high-dimensional
fixed effects.

{pstd}
The command combines the speed of the {cmd:creghdfe} HDFE absorption algorithm
(using conjugate gradient iteration) with comprehensive IV estimation. All heavy computation
is performed in optimized C code with OpenMP parallelization.

{pstd}
The syntax follows {cmd:ivreghdfe}: endogenous variables and their instruments
are specified in parentheses with an equals sign. The default estimator is 2SLS;
alternative estimators can be selected via options.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt absorb(varlist)} specifies categorical variables whose fixed effects are
absorbed from all variables (dependent, endogenous, exogenous, and instruments)
before IV estimation. This is required.

{dlgtab:Estimators}

{phang}
{opt liml} requests limited-information maximum likelihood estimation. LIML is
more robust to weak instruments than 2SLS but has the same asymptotic distribution
under standard conditions. The LIML k value (lambda) is stored in {cmd:e(lambda)}.

{phang}
{opt fuller(#)} requests Fuller's modified LIML estimator with adjustment parameter
{it:#}. Fuller(1) has minimum bias while Fuller(4) is approximately efficient. The
parameter must be positive. This estimator reduces the finite-sample bias of LIML.

{phang}
{opt kclass(#)} requests the k-class estimator with {it:k} = {it:#}. Special cases:
k=0 is OLS, k=1 is 2SLS, k=lambda is LIML. Values between 0 and 1 trade off between
OLS bias and 2SLS variance.

{phang}
{opt gmm2s} requests two-step efficient GMM estimation. This uses an optimal weighting
matrix computed from first-step residuals. More efficient than 2SLS under heteroskedasticity
when overidentified.

{phang}
{opt cue} requests continuously-updated GMM (CUE) estimation. CUE jointly estimates
coefficients and the optimal weighting matrix, providing better finite-sample
properties than two-step GMM.

{dlgtab:VCE/SE}

{phang}
{opt vce(vcetype)} specifies the variance-covariance estimation method:

{phang2}{opt unadjusted} or {opt ols} - standard VCE assuming homoskedasticity{p_end}
{phang2}{opt robust} - heteroskedasticity-robust standard errors (HC1){p_end}
{phang2}{opt cluster} {it:clustvar} - cluster-robust standard errors{p_end}

{pstd}
The VCE is computed using the proper sandwich formula for the chosen estimator
with the original endogenous regressors (not the first-stage fitted values).

{phang}
{opt small} requests small-sample corrections to be applied to variance estimates.
This includes adjusting for degrees of freedom lost due to absorbed fixed effects.

{dlgtab:HAC Standard Errors}

{phang}
{opt bw(#)} specifies the bandwidth for kernel-based heteroskedasticity and
autocorrelation consistent (HAC) standard errors. If not specified and a kernel
is requested, the Newey-West optimal bandwidth floor(4*(N/100)^(2/9)) is used.

{phang}
{opt kernel(string)} requests HAC standard errors with the specified kernel:

{phang2}{opt bartlett} or {opt nwest} - Bartlett (Newey-West) kernel{p_end}
{phang2}{opt parzen} - Parzen kernel{p_end}
{phang2}{opt quadraticspectral} or {opt qs} - Quadratic Spectral kernel{p_end}
{phang2}{opt truncated} or {opt trunc} - Truncated kernel{p_end}
{phang2}{opt tukey} or {opt thann} - Tukey-Hanning kernel{p_end}

{phang}
{opt dkraay(#)} requests Driscoll-Kraay standard errors with {it:#} lags, which
are robust to very general forms of cross-sectional and temporal dependence
in panel data. The Bartlett kernel is used by default.

{dlgtab:Estimation Settings}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the conjugate gradient
solver used in HDFE absorption. Default is 1e-8.

{phang}
{opt maxiter(#)} specifies the maximum number of CG iterations. Default is 500.

{phang}
{opt noconstant} suppresses the constant term. Note that with absorbed fixed
effects, the constant is typically absorbed anyway.

{dlgtab:Reporting}

{phang}
{opt first} reports first-stage regression statistics including F-statistics
for testing instrument strength for each endogenous variable.

{phang}
{opt rf} reports reduced-form regression estimates (regression of the dependent
variable directly on all instruments).

{phang}
{opt coviv} displays the covariance matrix of the IV estimators.

{phang}
{opt verbose} displays detailed progress information during computation.

{phang}
{opt timeit} displays a timing breakdown of computational phases.

{dlgtab:Save Results}

{phang}
{opt residuals(newvar)} saves the residuals from the IV regression to a new
variable named {it:newvar}.


{marker examples}{...}
{title:Examples}

{pstd}Basic 2SLS IV regression with one endogenous variable:{p_end}
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

{pstd}LIML estimation (more robust to weak instruments):{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) liml}{p_end}

{pstd}Fuller's modified LIML with alpha=1:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) fuller(1)}{p_end}

{pstd}Two-step GMM:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) gmm2s vce(robust)}{p_end}

{pstd}HAC standard errors with Bartlett kernel:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode) kernel(bartlett) bw(3)}{p_end}

{pstd}Driscoll-Kraay standard errors for panel data:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode) dkraay(4)}{p_end}

{pstd}With reduced-form and first-stage output:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) first rf}{p_end}


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
{synopt:{cmd:e(F)}}Wald F-statistic{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}

{p2col 5 20 24 2: First-stage scalars}{p_end}
{synopt:{cmd:e(F_first1)}}first-stage F-stat for 1st endogenous var{p_end}
{synopt:{cmd:e(F_first2)}}first-stage F-stat for 2nd endogenous var{p_end}
{synopt:{cmd:e(...)}}(additional F_first{it:k} for each endogenous var){p_end}

{p2col 5 20 24 2: Diagnostic test scalars}{p_end}
{synopt:{cmd:e(idstat)}}underidentification test statistic (Anderson or Kleibergen-Paap){p_end}
{synopt:{cmd:e(iddf)}}degrees of freedom for underidentification test{p_end}
{synopt:{cmd:e(idp)}}p-value for underidentification test{p_end}
{synopt:{cmd:e(cd_f)}}Cragg-Donald Wald F-statistic (weak identification){p_end}
{synopt:{cmd:e(kp_f)}}Kleibergen-Paap rk Wald F-statistic (if robust/cluster){p_end}
{synopt:{cmd:e(widstat)}}weak identification F-statistic (kp_f if robust, else cd_f){p_end}
{synopt:{cmd:e(sargan)}}Sargan or Hansen J statistic (overidentification){p_end}
{synopt:{cmd:e(sargan_df)}}degrees of freedom for overidentification test{p_end}
{synopt:{cmd:e(sargan_p)}}p-value for overidentification test{p_end}
{synopt:{cmd:e(endog_chi2)}}endogeneity test chi-squared statistic{p_end}
{synopt:{cmd:e(endog_f)}}endogeneity test F-statistic{p_end}
{synopt:{cmd:e(endog_df)}}degrees of freedom for endogeneity test{p_end}
{synopt:{cmd:e(endog_p)}}p-value for endogeneity test{p_end}

{p2col 5 20 24 2: LIML/Fuller scalars}{p_end}
{synopt:{cmd:e(lambda)}}LIML k-value (if liml or fuller){p_end}
{synopt:{cmd:e(fuller)}}Fuller parameter (if fuller){p_end}
{synopt:{cmd:e(kclass)}}k-class parameter (if kclass){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:civreghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(endogvars)}}endogenous variables{p_end}
{synopt:{cmd:e(exogvars)}}exogenous variables{p_end}
{synopt:{cmd:e(instruments)}}excluded instruments{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(title)}}estimation method title{p_end}
{synopt:{cmd:e(model)}}estimator: 2sls, liml, fuller, kclass, gmm2s, or cue{p_end}
{synopt:{cmd:e(vcetype)}}VCE type: Robust or Cluster{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustered){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}


{title:Diagnostic Tests}

{pstd}
{cmd:civreghdfe} automatically reports several diagnostic tests for IV estimation:

{dlgtab:Underidentification Test}

{pstd}
Tests whether the equation is identified (instruments are relevant). Uses the
Anderson canonical correlations LM statistic for unadjusted standard errors,
or the Kleibergen-Paap rk LM statistic for robust/cluster standard errors.
Rejection indicates instruments are relevant.

{dlgtab:Weak Identification Test}

{pstd}
The Cragg-Donald Wald F-statistic (or Kleibergen-Paap rk Wald F for robust/cluster)
tests whether instruments are weak. Compare to Stock-Yogo critical values shown
in output. Values above 10 generally indicate instruments are not weak.

{dlgtab:Overidentification Test}

{pstd}
The Sargan statistic (for unadjusted) or Hansen J statistic (for robust/cluster)
tests the joint null hypothesis that instruments are valid (uncorrelated with
the error term). Only available when overidentified. Rejection suggests some
instruments may be invalid.


{title:Methods and Formulas}

{pstd}
{cmd:civreghdfe} implements IV estimation with the following steps:

{phang2}1. Absorb fixed effects from all variables (y, endogenous X, exogenous X, instruments Z)
using the conjugate gradient (CG) solver with symmetric Kaczmarz iteration.{p_end}

{phang2}2. Compute the IV estimator on demeaned data:
{p_end}

{pmore}2SLS: {it:beta} = (X'P_Z X)^{-1} X'P_Z y where P_Z = Z(Z'Z)^{-1}Z'{p_end}

{pmore}LIML: Minimize eigenvalue lambda of (y,X)'M_Z(y,X) vs (y,X)'(y,X), then use k=lambda{p_end}

{pmore}Fuller: Use k = lambda - alpha/(N-K_iv) where alpha is the Fuller parameter{p_end}

{pmore}GMM: Use optimal weighting W = (Z'Omega Z)^{-1} where Omega is the residual covariance{p_end}

{phang2}3. Compute residuals using original X (not projected X) for proper inference.{p_end}

{phang2}4. Compute VCE using the appropriate sandwich estimator for the chosen method.{p_end}

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
sandwich estimator adapted for the chosen IV method.

{pstd}
Driscoll-Kraay standard errors are computed using the kernel-weighted
spatial correlation estimator, which is consistent under very general
forms of cross-sectional and temporal dependence.


{title:References}

{phang}
Anderson, T. W. and H. Rubin. 1949. Estimation of the parameters of a single
equation in a complete system of stochastic equations. {it:Annals of Mathematical Statistics} 20: 46-63.

{phang}
Cragg, J. G. and S. G. Donald. 1993. Testing identifiability and specification
in instrumental variable models. {it:Econometric Theory} 9: 222-240.

{phang}
Fuller, W. A. 1977. Some properties of a modification of the limited information
estimator. {it:Econometrica} 45: 939-953.

{phang}
Hansen, L. P. 1982. Large sample properties of generalized method of moments
estimators. {it:Econometrica} 50: 1029-1054.

{phang}
Kleibergen, F. and R. Paap. 2006. Generalized reduced rank tests using the
singular value decomposition. {it:Journal of Econometrics} 133: 97-126.

{phang}
Newey, W. K. and K. D. West. 1987. A simple, positive semi-definite, heteroskedasticity
and autocorrelation consistent covariance matrix. {it:Econometrica} 55: 703-708.

{phang}
Stock, J. H. and M. Yogo. 2005. Testing for weak instruments in linear IV regression.
In {it:Identification and Inference for Econometric Models: Essays in Honor of Thomas Rothenberg},
ed. D. W. K. Andrews and J. H. Stock, 80-108. Cambridge: Cambridge University Press.


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help ivreghdfe} (if installed), {help creghdfe}, {help ivreg2} (if installed),
{help reghdfe} (if installed)
{p_end}

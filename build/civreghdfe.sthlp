{smcl}
{* *! version 0.9.0 26Jan2026}{...}
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
{synopt:{opt vce(vcetype)}}variance-covariance estimation: {opt un:adjusted}, {opt r:obust}, {opt cl:uster} {it:clustvar} [{it:clustvar2}]{p_end}
{synopt:{opt r:obust}}heteroskedasticity-robust SEs (alias for vce(robust)){p_end}
{synopt:{opt cl:uster(varlist)}}cluster-robust SEs (alias for vce(cluster varlist)){p_end}
{synopt:{opt small}}use small-sample adjustments (df corrections){p_end}
{synopt:{opt dofminus(#)}}subtract {it:#} from residual degrees of freedom{p_end}
{synopt:{opt sdofminus(#)}}additional DOF adjustment for small-sample VCE correction{p_end}
{synopt:{opt nopartialsmall}}exclude partialled variables from K in small-sample adjustment{p_end}

{syntab:HAC Standard Errors}
{synopt:{opt bw(#)}}bandwidth for kernel-based HAC estimation{p_end}
{synopt:{opt kernel(string)}}kernel type for HAC: {opt bartlett}, {opt parzen}, {opt quadraticspectral}, {opt truncated}, {opt tukey}{p_end}
{synopt:{opt dkraay(#)}}Driscoll-Kraay SEs with {it:#} lags (for panel data){p_end}
{synopt:{opt center}}center score vectors before HAC outer product computation{p_end}

{syntab:Estimation Settings}
{synopt:{opt tol:erance(#)}}convergence tolerance for CG solver (default: 1e-8){p_end}
{synopt:{opt max:iter(#)}}maximum CG iterations (default: 500){p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt noc:onstant}}suppress constant term (absorbed with FE){p_end}

{syntab:Reporting}
{synopt:{opt first}}report first-stage regression statistics{p_end}
{synopt:{opt ff:irst}}report full first-stage statistics (partial R², F-stat){p_end}
{synopt:{opt rf}}report reduced-form estimates{p_end}
{synopt:{opt coviv}}display covariance matrix of IV estimators{p_end}
{synopt:{opt v:erbose}}display progress information{p_end}
{synopt:{opt timeit}}display timing breakdown{p_end}

{syntab:Display}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt nohead:er}}suppress output header{p_end}
{synopt:{opt nofoot:er}}suppress output footer (diagnostic tests){p_end}
{synopt:{opt noout:put}}suppress coefficient table{p_end}
{synopt:{opt ti:tle(string)}}custom title for output{p_end}
{synopt:{opt sub:title(string)}}subtitle displayed below title{p_end}
{synopt:{opt depn:ame(string)}}custom label for dependent variable in output{p_end}
{synopt:{opt noid}}suppress underidentification test display{p_end}
{synopt:{opt eform(string)}}report exponentiated coefficients with label {it:string}{p_end}
{synopt:{opt plus}}draw plus sign on coefficient table{p_end}
{synopt:{opt noomit:ted}}do not display omitted reference categories{p_end}
{synopt:{opt omit:ted}}display omitted reference categories{p_end}
{synopt:{opt vsquish}}suppress blank space in output{p_end}
{synopt:{opt noemptycells}}do not display empty cells for interactions{p_end}
{synopt:{opt baselev:els}}display base levels of factor variables{p_end}
{synopt:{opt allbaselev:els}}display all base levels{p_end}

{syntab:Diagnostic Tests}
{synopt:{opt orth:og(varlist)}}test orthogonality/exogeneity of specified excluded instruments (C-statistic){p_end}
{synopt:{opt endogt:est(varlist)}}test endogeneity of specified endogenous regressors{p_end}
{synopt:{opt red:undant(varlist)}}test redundancy of specified excluded instruments{p_end}

{syntab:Estimation}
{synopt:{opt part:ial(varlist)}}partial out specified exogenous regressors via FWL{p_end}
{synopt:{opt fwl(varlist)}}alias for {opt partial()} (Frisch-Waugh-Lovell){p_end}

{syntab:Save Results}
{synopt:{opt res:iduals(newvar)}}save residuals to new variable{p_end}
{synopt:{opt res:iduals2}}save residuals to auto-named variable _civreghdfe_resid{p_end}
{synopt:{opt savef:irst}}store estimation results using {cmd:estimates store}{p_end}
{synopt:{opt savefp:refix(string)}}prefix for stored results (default: _civreghdfe_){p_end}
{synopt:{opt saverf}}store reduced-form estimation results{p_end}
{synopt:{opt saverfp:refix(string)}}prefix for stored RF results (default: _civreghdfe_rf_){p_end}
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
{phang2}{opt cluster} {it:clustvar1} {it:clustvar2} - two-way cluster-robust standard errors{p_end}

{pstd}
Two-way clustering computes standard errors that are robust to arbitrary correlation
within each dimension using the Cameron-Gelbach-Miller (2011) formula:
V_twoway = V_clustvar1 + V_clustvar2 - V_intersection.
This is appropriate for panel data with both cross-sectional and time-series correlation.

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
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations including HDFE absorption and matrix computations. By default,
{cmd:civreghdfe} uses all available CPU cores as reported by OpenMP. Use this
option to limit parallelism, for example when running multiple jobs simultaneously.

{phang}
{opt noconstant} suppresses the constant term. Note that with absorbed fixed
effects, the constant is typically absorbed anyway.

{dlgtab:Reporting}

{phang}
{opt first} reports first-stage regression statistics including F-statistics
for testing instrument strength for each endogenous variable.

{phang}
{opt ffirst} reports full first-stage regression summary statistics in a
tabular format. For each endogenous variable, displays the partial R-squared
(the R² from regressing the endogenous variable on all instruments after
partialling out exogenous regressors) and the first-stage F-statistic with
its p-value. These statistics are also stored in {cmd:e()} for each
endogenous variable (e.g., {cmd:e(partial_r2_1)}, {cmd:e(F_first1)}).

{phang}
{opt rf} reports reduced-form regression estimates (regression of the dependent
variable directly on all instruments).

{phang}
{opt coviv} displays the covariance matrix of the IV estimators.

{phang}
{opt verbose} displays detailed progress information during computation.

{phang}
{opt timeit} displays a timing breakdown of computational phases.

{dlgtab:Display}

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals. The default is {cmd:level(95)} or as set by {helpb set level}.

{phang}
{opt noheader} suppresses the output header, which includes the estimation
method, sample size, R-squared, and other summary statistics.

{phang}
{opt nofooter} suppresses the output footer, which includes the diagnostic
tests (underidentification test, weak identification test, Sargan/Hansen test,
Stock-Yogo critical values, and instruments list).

{phang}
{opt nooutput} suppresses the coefficient table. The estimation is still
performed and results are stored in {cmd:e()}.

{phang}
{opt title(string)} replaces the default title "IV (2SLS) estimation" with
a custom title in the output.

{phang}
{opt depname(string)} specifies a custom label for the dependent variable
in the coefficient table output and stored results.

{phang}
{opt noid} suppresses the underidentification test display in the footer.
The test is still computed and stored in {cmd:e(idstat)}.

{dlgtab:Diagnostic Tests}

{phang}
{opt orthog(varlist)} requests the C-statistic (orthogonality test) for the specified
excluded instruments. The test evaluates whether the specified instruments satisfy
the exclusion restriction (i.e., are uncorrelated with the structural error).
The C-statistic is computed as C = J_full - J_restricted, where J_restricted is the
Sargan/Hansen statistic from a model excluding the tested instruments. Under H0 (that
the specified instruments are exogenous), C ~ chi-sq(df) where df = number of tested
instruments. A large C-statistic (small p-value) suggests the instruments may be
endogenous. Results are stored in {cmd:e(cstat)}, {cmd:e(cstat_df)}, and {cmd:e(cstat_p)}.

{phang}
{opt endogtest(varlist)} requests the endogeneity test for the specified endogenous
regressors. The test evaluates whether the specified variables can be treated as
exogenous. This uses a Durbin-Wu-Hausman style augmented regression test applied
to the subset of variables specified. Under H0 (that the specified regressors are
exogenous), the test statistic ~ chi-sq(df) where df = number of tested regressors.
A large test statistic (small p-value) supports treating these as endogenous.
Results are stored in {cmd:e(endogtest)}, {cmd:e(endogtest_df)}, and {cmd:e(endogtest_p)}.

{phang}
{opt redundant(varlist)} requests the instrument redundancy test for the specified
excluded instruments. The test evaluates whether the specified instruments contribute
to identification beyond the remaining instruments. Under H0 (that the instruments
are redundant), the LM statistic ~ chi-sq(df) where df = K_endog * (number tested).
A small test statistic (large p-value) suggests the instruments may be safely dropped.
Results are stored in {cmd:e(redund)}, {cmd:e(redund_df)}, and {cmd:e(redund_p)}.

{dlgtab:Estimation}

{phang}
{opt partial(varlist)} requests that the specified exogenous variables be partialled
out (removed) from the estimation using the Frisch-Waugh-Lovell (FWL) transformation.
All other variables (dependent, endogenous, remaining exogenous, and instruments)
are regressed on the partial variables and replaced with their residuals. The FWL
theorem guarantees that coefficients on the non-partialled variables are numerically
identical to those from the full regression. This option is useful when you want to
omit certain control variables from the displayed output while still controlling for
them. The partialled-out variables do not appear in the coefficient table.
Note: Variables specified in partial() must be included in the model as exogenous
regressors. The partialling is performed before fixed effect absorption.

{dlgtab:Save Results}

{phang}
{opt residuals(newvar)} saves the residuals from the IV regression to a new
variable named {it:newvar}.

{phang}
{opt savefirst} stores the current estimation results using {cmd:estimates store}
so they can be retrieved later with {cmd:estimates restore}. The results are
stored with a prefix (default: "_civreghdfe_") followed by "main".

{phang}
{opt savefprefix(string)} specifies the prefix to use when storing first-stage
results with {opt savefirst}. The default prefix is "_civreghdfe_".

{phang}
{opt saverf} stores the current estimation results (as reduced-form estimates)
using {cmd:estimates store}. The results are stored with a prefix (default:
"_civreghdfe_rf_") followed by "main".

{phang}
{opt saverfprefix(string)} specifies the prefix to use when storing reduced-form
results with {opt saverf}. The default prefix is "_civreghdfe_rf_".


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

{pstd}Test orthogonality of a specific instrument (C-statistic):{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue south) age, absorb(idcode) orthog(south)}{p_end}

{pstd}Test endogeneity of a specific regressor:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure hours = union wks_ue south) age, absorb(idcode) endogtest(hours)}{p_end}

{pstd}Test instrument redundancy:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue south) age, absorb(idcode) redundant(south)}{p_end}

{pstd}Partial out exogenous variables (FWL transformation):{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age ttl_exp, absorb(idcode) partial(ttl_exp)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:civreghdfe} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(df_a)}}absorbed degrees of freedom{p_end}
{synopt:{cmd:e(K)}}number of regressors{p_end}
{synopt:{cmd:e(K_endog)}}number of endogenous regressors{p_end}
{synopt:{cmd:e(K_exog)}}number of exogenous regressors{p_end}
{synopt:{cmd:e(K_iv)}}number of instruments (including exogenous){p_end}
{synopt:{cmd:e(G)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed FE dimensions (same as G){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if clustered); for two-way equals e(N_clust1){p_end}
{synopt:{cmd:e(N_clust1)}}number of first-dimension clusters (if two-way clustered){p_end}
{synopt:{cmd:e(N_clust2)}}number of second-dimension clusters (if two-way clustered){p_end}
{synopt:{cmd:e(F)}}Wald F-statistic{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: First-stage scalars}{p_end}
{synopt:{cmd:e(F_first1)}}first-stage F-stat for 1st endogenous var{p_end}
{synopt:{cmd:e(F_first2)}}first-stage F-stat for 2nd endogenous var{p_end}
{synopt:{cmd:e(...)}}(additional F_first{it:k} for each endogenous var){p_end}
{synopt:{cmd:e(partial_r2_1)}}partial R² for 1st endogenous var (if ffirst){p_end}
{synopt:{cmd:e(partial_r2_2)}}partial R² for 2nd endogenous var (if ffirst){p_end}

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

{p2col 5 20 24 2: Optional diagnostic test scalars}{p_end}
{synopt:{cmd:e(cstat)}}C-statistic for instrument orthogonality (if orthog() specified){p_end}
{synopt:{cmd:e(cstat_df)}}degrees of freedom for C-statistic{p_end}
{synopt:{cmd:e(cstat_p)}}p-value for C-statistic{p_end}
{synopt:{cmd:e(endogtest)}}endogeneity test for subset (if endogtest() specified){p_end}
{synopt:{cmd:e(endogtest_df)}}degrees of freedom for subset endogeneity test{p_end}
{synopt:{cmd:e(endogtest_p)}}p-value for subset endogeneity test{p_end}
{synopt:{cmd:e(redund)}}instrument redundancy test statistic (if redundant() specified){p_end}
{synopt:{cmd:e(redund_df)}}degrees of freedom for redundancy test{p_end}
{synopt:{cmd:e(redund_p)}}p-value for redundancy test{p_end}

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
{synopt:{cmd:e(vcetype)}}VCE type: Robust, Cluster, or Two-way Cluster{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustered); for two-way equals e(clustvar1){p_end}
{synopt:{cmd:e(clustvar1)}}first cluster variable (if two-way clustered){p_end}
{synopt:{cmd:e(clustvar2)}}second cluster variable (if two-way clustered){p_end}

{p2col 5 20 24 2: DOF adjustment scalars (if specified)}{p_end}
{synopt:{cmd:e(dofminus)}}user-specified DOF adjustment (if {cmd:dofminus()} specified){p_end}
{synopt:{cmd:e(sdofminus_opt)}}user-specified small-sample DOF adjustment (if {cmd:sdofminus()} specified){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


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
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command provides similar functionality to {cmd:ivreghdfe} by Sergio Correia.
The original {cmd:ivreghdfe} package, which combines {cmd:ivreg2} and {cmd:reghdfe},
is available at {browse "https://github.com/sergiocorreia/ivreghdfe"}.
We thank Sergio Correia for his groundbreaking work on high-dimensional fixed effects
estimation, as well as Mark Schaffer, Christopher Baum, and Steven Stillman for their
foundational work on {cmd:ivreg2}.

{pstd}
See also: Correia, S. 2016. "REGHDFE: Stata module to perform linear or instrumental-variable
regression absorbing any number of high-dimensional fixed effects." Statistical Software
Components, Boston College Department of Economics.


{title:Also see}

{psee}
Online: {browse "https://github.com/sergiocorreia/ivreghdfe":ivreghdfe} (if installed),
{browse "https://github.com/sergiocorreia/reghdfe":reghdfe} (if installed),
{help ivreg2} (if installed), {help creghdfe}
{p_end}

{smcl}
{* *! version 1.2.4 20261103}{...}
{viewerjumpto "Syntax" "civreghdfe##syntax"}{...}
{viewerjumpto "Description" "civreghdfe##description"}{...}
{viewerjumpto "Options" "civreghdfe##options"}{...}
{viewerjumpto "Prediction" "civreghdfe##prediction"}{...}
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
[{cmd:(}{it:endogvars} {cmd:=} {it:instruments}{cmd:)}]
[{it:exogvars}]
{ifin}
{weight}
[{cmd:,} {opt a:bsorb(varlist)} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt a:bsorb(varlist)}}categorical variables to absorb as fixed effects{p_end}

{syntab:Estimators}
{synopt:{opt liml}}use limited-information maximum likelihood (LIML){p_end}
{synopt:{opt fuller(#)}}use Fuller's modified LIML with parameter {it:#}{p_end}
{synopt:{opt kclass(#)}}use k-class estimator with {it:k} = {it:#}{p_end}
{synopt:{opt gmm2s}}use two-step efficient GMM{p_end}
{synopt:{opt cue}}use continuously-updated GMM (CUE){p_end}
{synopt:{opt b0(matrix)}}evaluate the CUE objective and statistics at the coefficients in {it:matrix}{p_end}

{syntab:VCE/SE}
{synopt:{opt vce(vcetype)}}variance-covariance estimation: {opt un:adjusted}, {opt r:obust}, {opt cl:uster} {it:clustvar} [{it:clustvar2}]{p_end}
{synopt:{opt r:obust}}heteroskedasticity-robust SEs (alias for vce(robust)){p_end}
{synopt:{opt cl:uster(varlist)}}cluster-robust SEs (alias for vce(cluster varlist)); may be combined with {opt robust}{p_end}
{synopt:{opt small}}use small-sample adjustments (df corrections){p_end}
{synopt:{opt dofminus(#)}}subtract {it:#} from residual degrees of freedom{p_end}
{synopt:{opt sdofminus(#)}}additional DOF adjustment for small-sample VCE correction{p_end}
{synopt:{opt nopartialsmall}}exclude partialled variables from K in small-sample adjustment{p_end}
{synopt:{opt dof:adjustments(list)}}absorbed degrees of freedom, as in {cmd:reghdfe}: {opt all} (default), {opt none}, {opt first:pair}, {opt pair:wise}, {opt cl:usters}{p_end}
{synopt:{opt keepsin:gletons}}keep singleton groups of the absorbed effects{p_end}

{syntab:HAC Standard Errors (data must be {help tsset})}
{synopt:{opt bw(#)}}bandwidth for kernel-based estimation: AC, or HAC with {opt robust}{p_end}
{synopt:{opt kernel(string)}}kernel (requires {opt bw()}): {opt bar:tlett} (default), {opt par:zen}, {opt qs}, {opt tru:ncated}, {opt thann}{p_end}
{synopt:{opt dkraay(#)}}Driscoll-Kraay SEs with bandwidth {it:#} (tsset panel data){p_end}
{synopt:{opt kiefer}}Kiefer SEs: within-panel autocorrelation (tsset panel data){p_end}

{syntab:Estimation Settings}
{synopt:{opt tol:erance(#)}}convergence tolerance for CG solver (default: 1e-8){p_end}
{synopt:{opt max:iter(#)}}maximum CG iterations (default: 16000){p_end}
{synopt:{opt iterate(#)}}synonym for {opt maxiter()}, as in {cmd:ivreghdfe}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt noc:onstant}}suppress constant term (absorbed with FE){p_end}
{synopt:{opt nocollin}}skip {cmd:ivreg2}'s collinearity checks of the variables{p_end}

{syntab:Reporting}
{synopt:{opt first}}report first-stage regression statistics{p_end}
{synopt:{opt ffir:st}}report full first-stage statistics (partial R², F-stat){p_end}
{synopt:{opt rf}}report reduced-form estimates{p_end}
{synopt:{opt coviv}}LIML and k-class covariance from the 2SLS form (X'P_Z X)^-1{p_end}
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}

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
{synopt:{opt coefl:egend}}display the coefficient legend instead of the table{p_end}
{synopt:{opt cformat(%fmt)}}format of coefficients, standard errors and confidence limits{p_end}
{synopt:{opt noemptycells}}do not display empty cells for interactions{p_end}
{synopt:{opt baselev:els}}display base levels of factor variables{p_end}
{synopt:{opt allbaselev:els}}display all base levels{p_end}

{syntab:Diagnostic Tests}
{synopt:{opt orth:og(varlist)}}test orthogonality/exogeneity of specified excluded instruments (C-statistic){p_end}
{synopt:{opt endog:test(varlist)}}test endogeneity of specified endogenous regressors{p_end}
{synopt:{opt red:undant(varlist)}}test redundancy of specified excluded instruments{p_end}

{syntab:Estimation}
{synopt:{opt part:ial(varlist)}}partial out specified exogenous regressors via FWL{p_end}
{synopt:{opt fwl(varlist)}}alias for {opt partial()} (Frisch-Waugh-Lovell){p_end}

{syntab:Save Results}
{synopt:{opt res:iduals(newvar)}}save residuals to new variable{p_end}
{synopt:{opt residuals2}}save residuals to auto-named variable _civreghdfe_resid{p_end}
{synopt:{opt savef:irst}}store the first-stage regressions using {cmd:estimates store}{p_end}
{synopt:{opt savefp:refix(string)}}prefix for stored results (default: _civreghdfe_){p_end}
{synopt:{opt saverf}}store reduced-form estimation results{p_end}
{synopt:{opt saverfp:refix(string)}}prefix for stored RF results (default: _civreghdfe_rf_){p_end}
{synopt:{opt noreturn}}skip all result storage (faster for scripted use){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{opt aweight}s, {opt fweight}s, and {opt pweight}s are allowed;
see {help weight}. {cmd:by} is allowed; see {help prefix}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:civreghdfe} estimates instrumental-variable models by 2SLS, LIML,
k-class, two-step GMM, or continuously updated GMM. C routines absorb
categorical fixed effects and compute the estimates and covariance matrix.


{pstd}
The syntax follows {cmd:ivreghdfe}: endogenous variables and their instruments
are specified in parentheses with an equals sign. Fixed effects are optional.
The default estimator is 2SLS; alternative estimators can be selected via
options. Without the parenthesized instrument specification, the command
estimates OLS. Exogenous regressors may appear before or after the
parentheses. Factor-variable and time-series operators are allowed in
regressor and instrument lists, and time-series operators on the dependent
variable. As in {cmd:ivreg2}, a term listed twice is used once, a term listed
both as exogenous and as endogenous is exogenous, and an excluded instrument
that is also a regressor is not an excluded instrument; these are listed in
{cmd:e(dups)}.

{pstd}
Collinear variables are handled as in {cmd:ivreg2}: a term collinear with earlier
terms is dropped ({cmd:e(collin)}), first among the endogenous regressors and then
over the exogenous regressors, excluded instruments and endogenous regressors in
that order; an endogenous regressor collinear with the instruments becomes an
included exogenous regressor ({cmd:e(ecollin)}), and the exogenous term it duplicates
is dropped. Warnings and the footer list them. {opt nocollin} skips these rules.

{pstd}
Differences from {cmd:ivreghdfe}: a regressor or instrument that the fixed
effects absorb is marked omitted ({cmd:o.}{it:varname}), where {cmd:ivreghdfe}
leaves its name unmarked with a zero coefficient and variance;
{cmd:civreghdfe} also omits variables whose remaining variation after removing
fixed effects falls below its numerical threshold. {cmd:ivreghdfe} may retain
these variables, producing different coefficients and degrees of freedom.
{opt sfirst} is not available.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt absorb(varlist)} specifies categorical variables whose fixed effects are
absorbed from all variables (dependent, endogenous, exogenous, and instruments)
before IV estimation. As in {cmd:reghdfe}, each term is a numeric or string
variable, an {cmd:i.}{it:var}, or an interaction {it:var1}{cmd:#}{it:var2}
absorbed as one effect; heterogeneous slopes and suboptions are not supported.
The cluster variables of {opt vce(cluster)} may likewise be string variables or
interactions.

{dlgtab:Estimators}

{phang}
{opt liml} requests limited-information maximum likelihood estimation. LIML is
computed from a generalized eigenvalue problem. The LIML k value (lambda) is stored in {cmd:e(lambda)}.

{phang}
{opt fuller(#)} requests Fuller's modified LIML estimator with adjustment parameter
{it:#}, k = lambda - {it:#}/(N - L). As in {cmd:ivreg2}, a negative {it:#} gives LIML.

{phang}
{opt kclass(#)} requests the k-class estimator with {it:k} = {it:#}, which must
not be negative. {cmd:kclass(1)} is 2SLS and k=lambda is LIML; {cmd:kclass(0)}
(the default) requests no k-class estimator.

{phang}
As in {cmd:ivreg2}, {opt fuller()} and {opt kclass()} take precedence over
{opt gmm2s}, and {opt liml} with {opt kclass()} is the k-class estimator with
that k (labelled LIML). {opt liml} with {opt gmm2s} or {opt cue}, {opt gmm2s}
with {opt cue}, a negative {opt kclass()}, {opt liml} with a negative
{opt fuller()}, {opt fuller()} with {opt kclass()} and {opt fuller()} or
{opt kclass()} with {opt cue} are errors (r(198)). These restrictions differ
from {cmd:ivreg2}: some combinations produce error 3200 there, while combining
{opt cue} with {opt fuller()} or {opt kclass()} runs CUE after a k-class fit.

{phang}
{opt gmm2s} requests two-step efficient GMM estimation. This uses an optimal weighting
matrix computed from first-step residuals.

{phang}
{opt cue} requests continuously-updated GMM (CUE) estimation. CUE recomputes the moment covariance at each candidate coefficient vector.

{phang}
{opt b0(matrix)} evaluates the model at the coefficients in {it:matrix}, as
{cmd:ivreg2} does: {cmd:e(b)} is {it:matrix}, {cmd:e(j)} the CUE objective (the J
statistic, also when the model is exactly identified), {cmd:e(V)} the efficient GMM
covariance at those coefficients, and the residual statistics follow from them
({cmd:e(model)} is {cmd:cue}). The columns of {it:matrix} are matched to the
regressors by name ({cmd:_cons} included without {opt absorb()}). {opt b0()} implies
{opt noid} and cannot be combined with another estimator.

{dlgtab:VCE/SE}

{phang}
{opt vce(vcetype)} specifies the variance-covariance estimation method:

{phang2}{opt unadjusted} or {opt ols} - standard VCE assuming homoskedasticity{p_end}
{phang2}{opt robust} - heteroskedasticity-robust standard errors{p_end}
{phang2}{opt cluster} {it:clustvar} - cluster-robust standard errors{p_end}
{phang2}{opt cluster} {it:clustvar1} {it:clustvar2} - two-way cluster-robust standard errors{p_end}

{pstd}
Two-way clustering computes standard errors that are robust to arbitrary correlation
within each dimension using the Cameron-Gelbach-Miller (2011) formula:
V_twoway = V_clustvar1 + V_clustvar2 - V_intersection.
This is appropriate for panel data with both cross-sectional and time-series correlation.

{phang}
{opt small} applies finite-sample covariance adjustments and reports t and F
inference. It is implied by {opt absorb()}. Without either option, coefficient
inference uses the normal distribution.

{phang}
{opt dofminus(#)} subtracts a nonnegative integer from the residual degrees
of freedom. {opt sdofminus(#)} supplies an additional nonnegative adjustment
for the small-sample covariance multiplier. {opt nopartialsmall} excludes
partialled-out regressors from that multiplier's regressor count.


{dlgtab:HAC Standard Errors}

{pstd}
The kernel-based estimators follow {cmd:ivreg2}. The data must be {help tsset}
(or {help xtset}) with a time variable; kernel lags are differences of the time
variable (in units of its delta) within the panels, so gaps and unbalanced panels
pair observations as {cmd:ivreg2}'s lag operators do, and the panel structure does
not depend on {opt absorb()}. fweights are not allowed with kernel-based estimators.

{phang}
{opt bw(#)} requests kernel-based standard errors and statistics with bandwidth
{it:#}, a positive real number as in {cmd:ivreg2} (lag windows sum the lags up to
{it:#}). Without {opt robust} the estimator is autocorrelation-consistent (AC) under
conditional homoskedasticity; with {opt robust} (or pweights) it is heteroskedasticity-
and autocorrelation-consistent (HAC). Combined with {cmd:cluster(}{it:timevar}{cmd:)}
it gives Driscoll-Kraay standard errors; combined with
{cmd:cluster(}{it:panelvar timevar}{cmd:)} two-way clustering on panel and time with
the kernel over time (panel clusters + Driscoll-Kraay - HAC). A kernel cannot be
combined with clustering on other variables. The bandwidth may not exceed the time
span of the estimation sample (in units of the tsset delta).

{phang}
{opt kernel(string)} sets the kernel and requires {opt bw()}:

{phang2}{opt bar:tlett} or {opt nwest} - Bartlett (Newey-West), the default{p_end}
{phang2}{opt par:zen} - Parzen{p_end}
{phang2}{opt qs}, {opt qua} or {opt quadraticspectral} - Quadratic Spectral (uses every lag up to the time span){p_end}
{phang2}{opt tru:ncated} or {opt trunc} - Truncated{p_end}
{phang2}{opt thann}, {opt tukey} or {opt tukey-hanning} - Tukey-Hanning{p_end}

{pstd}
{cmd:bw(auto)}, {opt center}, and the Tukey-Hamming, Daniell, and Tent kernels
are not supported.

{phang}
{opt dkraay(#)} requests Driscoll-Kraay standard errors, robust to heteroskedasticity
and to disturbances common to the panels that are autocorrelated: the clusters are
the time periods and the kernel ({opt bartlett} unless {opt kernel()} is specified)
with bandwidth {it:#} runs over the periods in calendar order. The data must be
tsset with panel and time variables. {opt dkraay(#)} is equivalent to
{cmd:cluster(}{it:timevar}{cmd:) bw(}{it:#}{cmd:)}.

{phang}
{opt kiefer} requests Kiefer (1980) standard errors, robust to arbitrary
within-panel autocorrelation but not to heteroskedasticity: the truncated kernel
with a bandwidth equal to the time span T of the sample, reported in {cmd:e(bw)}.
The data must be tsset with panel and time variables. {opt kiefer} cannot be
combined with {opt robust}, {opt cluster()}, {opt bw()} or {opt kernel()}; with
pweights the statistics are heteroskedasticity-robust as well, as in {cmd:ivreg2}.

{pstd}
The VCE, the J statistic, the C statistic, the endogeneity test and the efficient
GMM weighting matrix of {opt gmm2s} and {opt cue} all use the moment covariance of
the chosen VCE, including two-way clustering and the kernel estimators. As in
{cmd:ivreghdfe}, the underidentification, weak-identification and redundancy
statistics and the model re-estimated for {opt endogtest()} use the Bartlett
kernel whatever {opt kernel()} is, and with {opt kiefer} they are the iid
(pweights: heteroskedasticity-robust) statistics.

{dlgtab:Estimation Settings}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for the conjugate gradient
solver used in HDFE absorption. Default is 1e-8.

{phang}
{opt maxiter(#)} specifies the maximum number of CG iterations. Default is
16000, as in {cmd:reghdfe} and {cmd:ivreghdfe}; {opt iterate(#)} is a synonym.
Both {opt maxiter()} and {opt tolerance()} must be positive. Failure to remove
fixed effects within the iteration limit returns error 430.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt noconstant} suppresses the constant term. With {opt absorb()}, the constant is absorbed with the fixed effects.

{dlgtab:Reporting}

{phang}
{opt first} displays the first-stage regressions. Use {opt ffirst}
to display a formatted first-stage summary table.

{phang}
{opt ffirst} reports full first-stage regression summary statistics in a
tabular format. For each endogenous variable, displays the partial R-squared
(the R² from regressing the endogenous variable on all instruments after
partialling out exogenous regressors) and the first-stage F-statistic with
its p-value. As in {cmd:ivreg2}, the F-statistic tests the excluded instruments
under the model's covariance estimator (robust, cluster or kernel when one is
specified), and its denominator degrees of freedom are G-1 with G clusters.
These statistics are also stored in {cmd:e()} for each
endogenous variable (e.g., {cmd:e(partial_r2_1)}, {cmd:e(F_first1)}).

{phang}
{opt rf} reports reduced-form regression estimates (regression of the dependent
variable directly on all instruments).

{phang}
The first-stage and reduced-form regressions of {opt first}, {opt rf},
{opt savefirst} and {opt saverf} are run by {cmd:creghdfe}, which clusters in one
dimension only; with two-way clustering these options are an error.

{phang}
{opt coviv} computes the covariance of the LIML, Fuller and k-class estimators
with the 2SLS form: (X'P_Z X)^-1 replaces the k-class (X'(I-kM_Z)X)^-1 in the
classical covariance and in the robust, cluster and kernel sandwiches, as in
{cmd:ivreg2}. It has no effect on other estimators.

{phang}
{opt verbose} displays progress and elapsed times.

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
Stock-Yogo critical values, and instruments list) and the reduced form shown by
{opt rf}. The statistics are still stored in {cmd:e()}.

{phang}
{opt nooutput} suppresses the coefficient table. The estimation is still
performed and results are stored in {cmd:e()}.

{phang}
{opt title(string)} replaces the default title "IV (2SLS) estimation" with
a custom title in the output.

{phang}
{opt depname(string)} specifies a custom name for the dependent variable
stored in {cmd:e(depvar)}. As in {cmd:ivreg2} and {cmd:ivreghdfe}, the
coefficient table still shows the actual dependent variable.

{phang}
{opt noid} suppresses the underidentification test display in the footer.
The test is still computed and stored in {cmd:e(idstat)}.

{dlgtab:Diagnostic Tests}

{phang}
{opt orthog(varlist)} requests the C-statistic (orthogonality test) for the specified
excluded instruments. The test evaluates whether the specified instruments satisfy
the exclusion restriction (i.e., are uncorrelated with the structural error).
The C-statistic is computed as C = J_full - J_restricted, where J_restricted is the
Sargan/Hansen statistic of the efficient GMM fit without the tested instruments,
weighted by the corresponding block of the full model's moment covariance S (as
{cmd:ivreg2} does), so it is robust to the chosen VCE (robust, cluster, two-way,
kernel). Under H0 (that
the specified instruments are exogenous), C ~ chi-sq(df) where df = number of tested
instruments. A large C-statistic (small p-value) suggests the instruments may be
endogenous. Results are stored in {cmd:e(cstat)}, {cmd:e(cstat_df)}, and {cmd:e(cstat_p)}.

{phang}
{opt endogtest(varlist)} requests the endogeneity test for the specified endogenous
regressors. The test evaluates whether the specified variables can be treated as
exogenous. This uses a difference-in-Sargan (C-statistic) approach: the tested
endogenous variables are reclassified as exogenous, that model is re-estimated
(2SLS, LIML or two-step GMM as the main model; 2SLS for k-class and CUE, as in
{cmd:ivreg2}) with the moment covariance of the chosen VCE, and its C statistic
for the reclassified variables is computed. Under H0 (that the specified regressors are exogenous), the test
statistic ~ chi-sq(df) where df = number of tested regressors. A large test
statistic (small p-value) supports treating these as endogenous.
Results are stored in {cmd:e(endogtest)}, {cmd:e(endogtest_df)}, and {cmd:e(endogtest_p)}.

{pstd}
The varlists of {opt orthog()}, {opt endogtest()}, {opt redundant()} and
{opt partial()} may use factor-variable and time-series operators; they are
expanded and matched term by term against the model's (expanded) variables, as
{cmd:ivreg2} does. Repeated names are counted once and base levels are skipped.
{opt orthog()} and {opt redundant()} accept excluded instruments only. A listed
instrument or endogenous regressor that was dropped as collinear (or absorbed by
the fixed effects) cannot be tested, and the command exits with an error.

{phang}
{opt redundant(varlist)} requests the instrument redundancy test for the specified
excluded instruments. The test evaluates whether the specified instruments contribute
to identification beyond the remaining instruments: {cmd:ranktest}'s LM test (the
Kleibergen-Paap version under robust, cluster or kernel VCEs). Under H0 (that the instruments
are redundant), the LM statistic ~ chi-sq(df) where df = K_endog * (number tested).
A large p-value means the test does not reject redundancy.
Results are stored in {cmd:e(redund)}, {cmd:e(redund_df)}, and {cmd:e(redund_p)}.

{dlgtab:Estimation}

{phang}
{opt partial(varlist)} requests that the specified exogenous variables be
partialled out (removed) from the estimation using the Frisch-Waugh-Lovell
(FWL) transformation. All other variables (dependent, endogenous, remaining
exogenous, and instruments) are regressed on the partial variables and
replaced with their residuals. The variables listed in {opt partial()} remain
controls in the model, but their coefficients are omitted from the coefficient
table. Variables specified in {opt partial()} must be included in the model as
exogenous regressors. The partialling is performed after fixed effect
absorption, on the demeaned data, and uses the estimation weights. A partial()
variable that is collinear with earlier ones or absorbed by the fixed effects
is dropped with a note. Without {opt absorb()}, the constant is partialled out
as well and counted in {cmd:e(partial_ct)}, as in {cmd:ivreg2}; {cmd:e(b)}
then has no {cmd:_cons}.

{dlgtab:Save Results}

{phang}
{opt residuals(newvar)} saves the residuals from the IV regression to a new
variable named {it:newvar}.

{phang}
{opt savefirst} stores each first-stage regression, as {cmd:ivreg2} does, with
{cmd:estimates store} under the prefix (default: "_civreghdfe_") followed by the
name of the endogenous regressor, so they can be retrieved later with
{cmd:estimates restore}. The constant is partialled out with the fixed effects.
The IV estimation results are also stored, under the prefix followed by "main".

{phang}
{opt savefprefix(string)} specifies the prefix to use when storing first-stage
results with {opt savefirst}. The default prefix is "_civreghdfe_".

{phang}
{opt saverf} stores the reduced-form regression of the dependent variable on
the instruments, as {cmd:ivreg2} does, with {cmd:estimates store} under the
prefix (default: "_civreghdfe_rf_") followed by the name of the dependent
variable. The IV estimation results are also stored, under the prefix followed
by "main".

{phang}
{opt saverfprefix(string)} specifies the prefix to use when storing reduced-form
results with {opt saverf}. The default prefix is "_civreghdfe_rf_".

{phang}
{opt noreturn} skips storing estimation results in {cmd:e()}. Use it when only
generated variables, such as residuals, are needed.


{marker prediction}{...}
{title:Prediction}

{pstd}
{cmd:predict newvar} defaults to {opt xb}, the structural linear prediction
excluding absorbed effects. {opt stdp} gives its standard error.
{opt residuals} returns estimation-sample residuals; with absorbed effects,
specify {opt residuals(newvar)} during estimation first. The option
{opt residuals2} saves these as {cmd:_civreghdfe_resid}.

{pstd}
{cmd:predict newvar, xbd} returns y minus the saved residual, including absorbed
effects. It requires saved residuals and is missing outside {cmd:e(sample)}.
{opt d} returns the sum of the absorbed effects, y - xb - the saved residual, and
{opt dresiduals} y - xb, as {cmd:reghdfe}'s {cmd:predict} does; both require
{opt residuals()} during estimation.

{pstd}
After {opt partial()}, {cmd:e(b)} omits the partialled-out terms
({cmd:e(partial)}, {cmd:e(partial_ct)}). Without {opt absorb()}, as with
{cmd:ivreg2}, only {opt residuals} can be predicted: y - Xb net of its
projection on the partialled-out terms, which equals the full model's
residual. With {opt absorb()}, {opt xb} and {opt stdp} refer to the posted
coefficients, as with {cmd:ivreghdfe}.

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

{pstd}LIML estimation:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) liml}{p_end}

{pstd}Fuller's modified LIML with alpha=1:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) fuller(1)}{p_end}

{pstd}Two-step GMM:{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) gmm2s vce(robust)}{p_end}

{pstd}HAC standard errors with Bartlett kernel:{p_end}
{phang2}{cmd:. xtset idcode year}{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode) kernel(bartlett) bw(3) robust}{p_end}

{pstd}AC standard errors (autocorrelation-consistent, homoskedastic):{p_end}
{phang2}{cmd:. civreghdfe ln_wage (tenure = union) age, absorb(idcode) bw(3)}{p_end}

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
{synopt:{cmd:e(df_m)}}model degrees of freedom, as {cmd:ivreg2}: rankxx - cons + sdofminus - partialcons{p_end}
{synopt:{cmd:e(rankxx)}, {cmd:e(rankzz)}}regressors and instruments kept (ranks){p_end}
{synopt:{cmd:e(rankS)}, {cmd:e(rankV)}}rank of the moment covariance S and nonzero variances in e(V){p_end}
{synopt:{cmd:e(sdofminus)}, {cmd:e(partialcons)}, {cmd:e(cons)}}{cmd:ivreg2}'s small-sample count of partialled-out terms, partialled-out constant (absorbed effects) and constant indicator{p_end}
{synopt:{cmd:e(ll)}}log likelihood, as {cmd:ivreg2} ({cmd:estat ic} uses it){p_end}
{synopt:{cmd:e(num_singletons)}, {cmd:e(N_full)}}singletons dropped and observations before they are dropped{p_end}
{synopt:{cmd:e(nocollin)}}1 if {opt nocollin} was specified{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom (with {opt small} or {opt absorb()}){p_end}
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
{synopt:{cmd:e(bw)}}bandwidth of a kernel-based VCE (Kiefer: the time span T){p_end}
{synopt:{cmd:e(F)}}Wald F-statistic{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: First-stage scalars}{p_end}
{synopt:{cmd:e(F_first1)}}first-stage F-stat for 1st endogenous var (missing for an endogenous variable dropped as collinear){p_end}
{synopt:{cmd:e(F_first2)}}first-stage F-stat for 2nd endogenous var{p_end}
{synopt:{cmd:e(...)}}(additional F_first{it:k} for each endogenous var){p_end}
{synopt:{cmd:e(partial_r2_1)}}partial R² for 1st endogenous var (if ffirst){p_end}
{synopt:{cmd:e(partial_r2_2)}}partial R² for 2nd endogenous var (if ffirst){p_end}

{p2col 5 20 24 2: Diagnostic test scalars}{p_end}
{synopt:{cmd:e(idstat)}}underidentification test statistic (Anderson or Kleibergen-Paap){p_end}
{synopt:{cmd:e(iddf)}}degrees of freedom for underidentification test{p_end}
{synopt:{cmd:e(idp)}}p-value for underidentification test{p_end}
{synopt:{cmd:e(cd_f)}}Cragg-Donald Wald F-statistic (weak identification){p_end}
{synopt:{cmd:e(kp_f)}}Kleibergen-Paap rk Wald F-statistic (if robust, cluster or kernel VCE){p_end}
{synopt:{cmd:e(widstat)}}weak identification F-statistic (kp_f if robust, cluster or kernel VCE, else cd_f){p_end}
{synopt:{cmd:e(sargan)}}Sargan or Hansen J statistic (overidentification; missing if S is rank deficient){p_end}
{synopt:{cmd:e(sargan_df)}}degrees of freedom for overidentification test{p_end}
{synopt:{cmd:e(sargan_p)}}p-value for overidentification test{p_end}
{synopt:{cmd:e(j)}, {cmd:e(jdf)}, {cmd:e(jp)}}the same, under {cmd:ivreg2}'s names{p_end}
{synopt:{cmd:e(cdf)}}the Cragg-Donald F, under {cmd:ivreg2}'s name{p_end}
{synopt:{cmd:e(Fp)}}p-value of {cmd:e(F)}{p_end}
{synopt:{cmd:e(r2c)}}centered R-squared (as {cmd:e(r2)}){p_end}
{synopt:{cmd:e(endog_ct)}, {cmd:e(exexog_ct)}, {cmd:e(inexog_ct)}}numbers of endogenous regressors, excluded instruments and included exogenous regressors{p_end}
{synopt:{cmd:e(endog_chi2)}}regression-based Wu-Hausman endogeneity test under iid errors: the Wald statistic of the first-stage residuals added to the (weighted) model, {cmd:e(endog_df)} times {cmd:e(endog_f)}{p_end}
{synopt:{cmd:e(endog_f)}}its F form, the F of {cmd:regress} {it:y X v}{cmd:; test} {it:v} for the first-stage residuals {it:v}{p_end}
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
{synopt:{cmd:e(kclass)}}k of the k-class estimator (if liml, fuller or kclass){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:civreghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(endogvars)}}endogenous variables{p_end}
{synopt:{cmd:e(exogvars)}}exogenous variables{p_end}
{synopt:{cmd:e(instruments)}}excluded instruments{p_end}
{synopt:{cmd:e(instd)}, {cmd:e(inexog)}, {cmd:e(exexog)}, {cmd:e(insts)}}instrumented, included exogenous, excluded and all instruments (expanded terms), as {cmd:ivreg2} names them{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(dups)}, {cmd:e(collin)}, {cmd:e(ecollin)}}duplicate terms, terms dropped as collinear, endogenous regressors reclassified as exogenous{p_end}
{synopt:{cmd:e(title)}}estimation method title{p_end}
{synopt:{cmd:e(model)}}estimator: 2sls, liml, fuller, kclass, gmm2s, or cue{p_end}
{synopt:{cmd:e(vcetype)}}VCE type: Robust, Cluster, or Two-way Cluster (Robust for Driscoll-Kraay; empty for AC and Kiefer){p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if clustered); contains both variable names for two-way clustering{p_end}
{synopt:{cmd:e(clustvar1)}}first cluster variable (if two-way clustered){p_end}
{synopt:{cmd:e(clustvar2)}}second cluster variable (if two-way clustered){p_end}
{synopt:{cmd:e(kernel)}}kernel of a kernel-based VCE: Bartlett, Parzen, Quadratic Spectral, Truncated or Tukey-Hanning{p_end}
{synopt:{cmd:e(tvar)}}tsset time variable (kernel-based VCE){p_end}
{synopt:{cmd:e(ivar)}}tsset panel variable (kernel-based VCE with panel data){p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression (if weighted){p_end}
{synopt:{cmd:e(marginsok)}}predictions allowed by {cmd:margins}{p_end}
{synopt:{cmd:e(marginsnotok)}}predictions disallowed by {cmd:margins} (with endogenous regressors){p_end}

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
or the Kleibergen-Paap rk LM statistic (for any number of endogenous regressors)
for robust, cluster, two-way, AC/HAC and Driscoll-Kraay standard errors, computed
as {cmd:ranktest} does. Rejection indicates instruments are relevant.

{dlgtab:Weak Identification Test}

{pstd}
The Cragg-Donald Wald F-statistic (or Kleibergen-Paap rk Wald F for robust,
cluster and kernel VCEs) tests whether instruments are weak. The output displays Stock-Yogo critical values for the Cragg-Donald statistic.
Those values assume independent, homoskedastic errors and do not provide
critical values for the robust Kleibergen-Paap statistic.

{dlgtab:Overidentification Test}

{pstd}
The Sargan statistic (for unadjusted, AC and Kiefer) or Hansen J statistic (for
robust/cluster) tests the joint null hypothesis that instruments are valid
(uncorrelated with the error term), with the moment covariance S of the chosen VCE
(for an inefficient estimator, the J of efficient GMM with that S, as in
{cmd:ivreg2}). Only available when overidentified. Rejection suggests some
instruments may be invalid.

{pstd}
When S is not of full rank (Stata's {cmd:invsym()} rank; typically too few clusters),
the J, C and endogeneity statistics are not reported (missing), as in {cmd:ivreg2},
and {opt gmm2s} and {opt cue} exit with error 506 because the optimal weighting
matrix is not unique.


{title:Methods and Formulas}

{pstd}
{cmd:civreghdfe} implements IV estimation with the following steps:

{phang2}1. Absorb fixed effects from all variables (y, endogenous X, exogenous X, instruments Z)
using the conjugate gradient (CG) solver with symmetric Kaczmarz iteration.
As in {cmd:ivreghdfe}, which partials the fixed effects out before it estimates,
the ranks are those of Stata's {cmd:invsym()} of the partialled (weighted)
cross-products: a regressor or instrument is dropped when the fixed effects absorb it
(its partialled sum of squares is rounding noise, at most 1e-20 times its sum of
squares) or when {cmd:invsym()} zeroes it; a regressor with little within variation
is kept. A dropped column does not count in K or the degrees of freedom; an absorbed
exogenous regressor or excluded instrument still counts among the instruments of the
Cragg-Donald and Kleibergen-Paap Wald F statistics and the identification dof, as
{cmd:ivreghdfe} counts its variable lists. An absorbed endogenous regressor has no
first-stage variation, so the identification statistics are 0. The footer reports the
instruments dropped this way.{p_end}

{phang2}2. Compute the IV estimator on demeaned data:
{p_end}

{pmore}2SLS: {it:beta} = (X'P_Z X)^{-1} X'P_Z y where P_Z = Z(Z'Z)^{-1}Z'{p_end}

{pmore}LIML: Minimize eigenvalue lambda of (y,X)'M_Z(y,X) vs (y,X)'(y,X), then use k=lambda{p_end}

{pmore}Fuller: Use k = lambda - alpha/(N-K_iv) where alpha is the Fuller parameter{p_end}

{pmore}GMM: Use optimal weighting W = S^{-1}, where S is the covariance of the moment conditions of the chosen VCE (robust, cluster, two-way, AC/HAC, Driscoll-Kraay, Kiefer) at the first-step (2SLS) residuals; CUE evaluates S at each candidate estimate{p_end}

{phang2}3. Compute residuals using the regressors rather than their
instrument-fitted values, after removing any absorbed fixed effects.{p_end}

{phang2}4. Compute VCE using the appropriate sandwich estimator for the chosen method.{p_end}

{pstd}
First-stage F-statistics are computed from the partial R-squared of the excluded
instruments in the first-stage regression.


{title:References}

{phang}
Anderson, T. W. and H. Rubin. 1949. Estimation of the parameters of a single
equation in a complete system of stochastic equations. {it:Annals of Mathematical Statistics} 20: 46-63.

{phang}
Cragg, J. G. and S. G. Donald. 1993. Testing identifiability and specification
in instrumental variable models. {it:Econometric Theory} 9: 222-240.

{phang}
Driscoll, J. C. and A. C. Kraay. 1998. Consistent covariance matrix estimation with
spatially dependent panel data. {it:Review of Economics and Statistics} 80: 549-560.

{phang}
Fuller, W. A. 1977. Some properties of a modification of the limited information
estimator. {it:Econometrica} 45: 939-953.

{phang}
Hansen, L. P. 1982. Large sample properties of generalized method of moments
estimators. {it:Econometrica} 50: 1029-1054.

{phang}
Kiefer, N. M. 1980. Estimation of fixed effect models for time series of cross-sections
with arbitrary intertemporal covariance. {it:Journal of Econometrics} 14: 195-202.

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
{cmd:civreghdfe} estimates instrumental-variable regressions with fixed
effects, similar to Sergio Correia's {cmd:ivreghdfe}. The original
{cmd:ivreghdfe} package, which combines {cmd:ivreg2} and {cmd:reghdfe}, is
available at {browse "https://github.com/sergiocorreia/ivreghdfe"}.
{cmd:ivreg2} is by Christopher Baum, Mark Schaffer, and Steven Stillman.

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

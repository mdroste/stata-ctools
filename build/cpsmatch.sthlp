{smcl}
{* *! version 1.2.6 20261111}{...}
{viewerjumpto "Syntax" "cpsmatch##syntax"}{...}
{viewerjumpto "Description" "cpsmatch##description"}{...}
{viewerjumpto "Options" "cpsmatch##options"}{...}
{viewerjumpto "Remarks" "cpsmatch##remarks"}{...}
{viewerjumpto "Examples" "cpsmatch##examples"}{...}
{viewerjumpto "Stored results" "cpsmatch##results"}{...}
{title:Title}

{phang}
{bf:cpsmatch} {hline 2} C-accelerated propensity score and Mahalanobis matching for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cpsmatch}
{it:treatvar}
[{it:varlist}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Propensity Score}
{synopt:{opth ps:core(varname)}}use existing propensity score variable{p_end}
{synopt:{opt logit}}estimate propensity score using logit{p_end}
{synopt:{opt probit}}estimate propensity score using probit (default){p_end}
{synopt:{opt index}}match on the linear index rather than the probability{p_end}
{synopt:{opt odds}}match on the log odds of the propensity score{p_end}
{synopt:{opt qui:etly}}do not display the score model{p_end}
{synopt:{opt warn:ings}}warn about controls with equal propensity scores{p_end}

{syntab:Matching Method}
{synopt:{opth n:eighbor(#)}}number of nearest neighbors (default: 1){p_end}
{synopt:{opth cal:iper(#)}}caliper width for matching{p_end}
{synopt:{opt radius}}radius matching{p_end}
{synopt:{opt kernel}}kernel matching{p_end}
{synopt:{opt llr}}local linear regression matching{p_end}
{synopt:{opth kerneltype(string)}}kernel function: epan, normal, biweight, uniform, tricube{p_end}
{synopt:{opth bw:idth(#)}}bandwidth for kernel and llr matching (default: 0.06){p_end}
{synopt:{opt spline}}spline-smoothing matching (needs the {cmd:spline} package){p_end}
{synopt:{opth nk:nots(#)}}number of interior knots for {opt spline}{p_end}
{synopt:{opth mahal:anobis(varlist)}}match on the Mahalanobis metric of {it:varlist}{p_end}
{synopt:{opt w(matrix)}}weighting matrix of the quadratic metric{p_end}

{syntab:Matching Options}
{synopt:{opt com:mon}}impose common support{p_end}
{synopt:{opt trim(#)}}impose common support by trimming # percent{p_end}
{synopt:{opt norep:lacement}}match without replacement{p_end}
{synopt:{opt ties}}include all tied matches{p_end}
{synopt:{opt desc:ending}}match in descending propensity score order{p_end}

{syntab:Outcome Analysis}
{synopt:{opth out:come(varlist)}}outcome variables{p_end}
{synopt:{opt ate}}also report the ATU and the ATE{p_end}
{synopt:{opt ai(#)}}Abadie-Imbens standard errors with # neighbors{p_end}
{synopt:{opt pop:ulation}}with {opt ai()}, the population variance{p_end}
{synopt:{opt altv:ariance}}with {opt ai()}, the conditional variance of Abadie et al. (2004){p_end}

{syntab:Performance}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cpsmatch} implements the matching estimators of {cmd:psmatch2} (version
4.0.12): nearest-neighbor, radius, kernel, local linear regression and
spline matching on a propensity score, and nearest-neighbor, kernel and local
linear regression matching on the Mahalanobis metric of covariates. With
{opt outcome()}, it reports the average treatment effect on the treated
(ATT), and with {opt ate} also on the untreated (ATU) and overall (ATE).
Observations with {cmd:_treated} equal to 1 are matched to those with
{cmd:_treated} equal to 0; see {help cpsmatch##remarks:Remarks} for other
treatment codes. Supply covariates for score estimation or specify
{opt pscore()}, not both. Covariates may include factor variables
({help fvvarlist}).

{pstd}
Matching on the propensity score (nearest neighbor, radius and kernel) runs
in a C plugin; Mahalanobis-metric matching and local linear regression
matching with a kernel other than the Epanechnikov kernel run in Mata, as
in {cmd:psmatch2}, and take time proportional to the product of the numbers
of treated and controls.

{pstd}
{cmd:cpsmatch} creates the variables {cmd:psmatch2} creates; see
{help cpsmatch##results:Stored results}.


{marker options}{...}
{title:Options}

{dlgtab:Propensity Score}

{phang}
{opth pscore(varname)} specifies an existing variable containing propensity
scores. If specified, {cmd:cpsmatch} uses this variable directly instead of
estimating the propensity score.

{phang}
{opt logit} specifies that the propensity score should be estimated using
logistic regression.

{phang}
{opt probit} specifies that the propensity score should be estimated using
probit regression. This is the default.

{phang}
{opt index} matches on the score model's linear index (its {cmd:predict, xb})
instead of the probability.

{phang}
{opt odds} matches on the log odds ln(p/(1-p)) of the score, stored in
{cmd:_pscore}. Common support ({opt common}, {opt trim()}) uses the score
itself. Observations without log odds (a score of 0 or 1) leave the sample.

{phang}
{opt quietly} does not display the output of the score model.

{phang}
{opt warnings} checks for supported controls with equal propensity scores
(with {opt ate}, also treated) and, with nearest-neighbor matching, warns that
the sort order of the data could affect the results.

{dlgtab:Matching Method}

{phang}
{opth neighbor(#)} specifies the number of nearest neighbors to match to each
treated observation. The default is 1. With {opt ties}, additional matches
may be included if there are ties at the boundary.

{phang}
{opth caliper(#)} specifies the maximum distance for a valid match:
observations without a match within the caliper are off support. On the
propensity score a match needs a distance strictly below the caliper; on the
Mahalanobis metric, a squared distance at most the caliper. Kernel and llr
matching ignore {opt caliper()}, as in {cmd:psmatch2}.

{phang}
{opt radius} specifies radius matching: every control with a score distance
at most {opt caliper()} is a match, with equal weights. When {opt caliper()} is
not specified or is 0, every control is a match, as in {cmd:psmatch2}, whose
missing bandwidth then spreads each match over every observation of the data:
as there, each on-support observation's {cmd:_weight} gains the number of
matched observations divided by {cmd:_N} (so the ATT is the plain difference
in means, and its standard error is {cmd:psmatch2}'s). A negative
{opt caliper()} matches nothing. With
both {opt radius} and {opt kernel}, radius matching is performed, as in
{cmd:psmatch2}.

{phang}
{opt kernel} specifies kernel matching: each control within the bandwidth
receives a weight from a kernel function of its distance to the treated
observation, normalized to sum to one. Compact-support kernels match nothing
when {opt bwidth()} is 0 or negative; the normal kernel has unbounded support
and uses the absolute value of {opt bwidth()}, as in {cmd:psmatch2}.

{phang}
{opt llr} specifies local linear regression matching: the matched outcome is
the intercept of a kernel-weighted regression of the controls' outcomes on
their distance. With the Epanechnikov kernel (the default) on the propensity
score, as {cmd:psmatch2} does, the controls' outcomes are smoothed with
{helpb lpoly} ({cmd:degree(1)}, bandwidth {opt bwidth()}) and stored in
{cmd:_s_}{it:outcome}, and each treated observation is matched to its nearest
neighbor(s)' smoothed values (so {cmd:_id}, {cmd:_n1} and {cmd:_nn} are
created). With another kernel or on the Mahalanobis metric, each treated
observation is matched by its own weighted regression; observations where it
cannot be estimated are off support. {opt llr} needs {opt outcome()} with the
Epanechnikov kernel on the score, and cannot be combined with {opt kernel}.

{phang}
{opth kerneltype(string)} specifies the kernel function for kernel and llr matching:

{p 12 16 2}{bf:epan} - Epanechnikov kernel (default){p_end}
{p 12 16 2}{bf:normal} - Gaussian/normal kernel{p_end}
{p 12 16 2}{bf:biweight} - Biweight kernel{p_end}
{p 12 16 2}{bf:uniform} - Uniform kernel{p_end}
{p 12 16 2}{bf:tricube} - Tricube kernel{p_end}

{phang}
{opth bwidth(#)} specifies the bandwidth for kernel and llr matching. The
default is 0.06, also for llr with the Epanechnikov kernel (as
{cmd:psmatch2} passes it to {cmd:lpoly}).

{phang}
{opt spline} performs spline-smoothing matching as {cmd:psmatch2} does: the
{cmd:spline} command (a user-written package) fits a natural cubic spline of
each outcome on the score among the supported controls, stored in
{cmd:_s_}{it:outcome}, and treated observations are matched to their nearest
neighbors' smoothed values. Without the package, {cmd:cpsmatch} exits with
error 198, as {cmd:psmatch2} does. The ATT standard error is missing.

{phang}
{opth nknots(#)} sets the number of interior knots for {opt spline}; the
default is the fourth root of the number of supported controls.

{phang}
{opth mahalanobis(varlist)} matches on the Mahalanobis distance of
{it:varlist} (factor variables allowed; collinear terms are dropped) instead
of on the propensity score, with the inverse of the pooled within-group
covariance matrix of the treated and controls. A propensity score, if
estimated or given, is stored but not used for matching unless it is in
{it:varlist}. {opt ties} is ignored, with a note, and {opt noreplacement}
is not allowed.

{phang}
{opt w(matrix)} specifies the weighting matrix W of the quadratic metric
(x-y)'W(x-y), in place of the inverse covariance matrix; its rows and
columns follow {opt mahalanobis()}.

{dlgtab:Matching Options}

{phang}
{opt common} imposes the common support restriction, excluding treated
observations with propensity scores outside the range of control propensity
scores (with {opt ate}, also controls outside the range of the treated).
Treated observations without an eligible match also receive
{cmd:_support=0}. It needs a propensity score.

{phang}
{opt trim(#)} imposes common support by dropping the treated observations at
which the density of the controls' propensity score (a {helpb kdensity}
estimate at each score) is zero or among its lowest # percent (with
{opt ate}, the controls likewise against the treated's density). # is
between 0 and 100; {cmd:trim(0)} drops only zero densities.

{phang}
{opt noreplacement} specifies matching without replacement. Each control can
only be matched to one treated observation. As in {cmd:psmatch2}, only
nearest-neighbor matching on the propensity score with {cmd:neighbor(1)}
is supported (with {opt radius} the option is ignored). Treated observations
are matched sequentially in propensity-score order (see {opt descending}). For
each one, the search starts at the next unused control after the previous
match and continues until the score distance increases. After reaching the
last control, the search reverses direction. Treated observations remaining
when no further control can be selected receive {cmd:_support=0}. With
replacement or without, the nearest control is the first of equally distant
controls in score order (the lower score on equal distance), as in
{cmd:psmatch2}.

{phang}
{opt ties} specifies that when there are tied matches at the boundary
distance, all tied observations should be included as matches.

{phang}
{opt descending} specifies that matching should proceed in descending
order of propensity score. As in {cmd:psmatch2}, it requires
{opt noreplacement}.

{dlgtab:Outcome Analysis}

{phang}
{opth outcome(varlist)} specifies the outcome variables. {cmd:cpsmatch}
computes and displays the ATT of each. As in {cmd:psmatch2}, the outcomes
share one match over the observations where all of them are nonmissing;
{cmd:r(att_}{it:var}{cmd:)} and {cmd:r(seatt_}{it:var}{cmd:)} hold each
outcome's results, and {cmd:r(att)} and {cmd:r(seatt)} the last one's.

{phang}
{opt ate} also matches each control to the treated and reports the average
treatment effect on the untreated (ATU) and the average treatment effect
(ATE), the ATT and the ATU weighted by the numbers of supported treated and
controls. Controls then receive matched outcomes, {cmd:_n#}, {cmd:_nn} and
{cmd:_pdif}; the treated's {cmd:_weight} counts their use as matches; and
common support and calipers apply to the controls as well. Without
{opt ai()}, the ATU and ATE standard errors are missing.

{phang}
{opt ai(#)} computes the heteroskedasticity-consistent standard errors of
Abadie and Imbens (2006, formula 14) with nearest-neighbor matching: each
observation is matched to its # nearest neighbors in its own arm, stored in
{cmd:_self_}{it:outcome}, to estimate the conditional outcome variance.
{cmd:psmatch2} documents them for Mahalanobis-metric matching; it also
computes them for nearest-neighbor matching on the score, and so does
{cmd:cpsmatch}.

{phang}
{opt population} with {opt ai()} reports the standard errors of the
population rather than the sample treatment effects.

{phang}
{opt altvariance} with {opt ai()} estimates the conditional variance as in
Abadie et al. (2004, p. 303). ({cmd:psmatch2} fails with this option when
there are several outcomes; {cmd:cpsmatch} does not.)

{dlgtab:Performance}

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
The matching sample excludes missing treatment, covariates, propensity scores,
Mahalanobis variables and, when specified, outcomes. Estimated scores also
restrict matching to the score model's estimation sample.

{pstd}
As in {cmd:psmatch2}, the treatment needs no particular coding. {cmd:_treated}
is the treatment stored as a byte: observations with {cmd:_treated} equal to 1
are matched to those equal to 0. Observations with another code stay in the
sample, on support, with {cmd:_weight} equal to {cmd:_treated}, and take no
part in matching. With no treated or no control observations the command runs
and the ATT is missing. {cmd:psmatch2} assumes every other code sorts after the
controls; a code below zero enters its matching as a control or a treated
observation, which {cmd:cpsmatch} does not do.

{pstd}
Each run replaces existing {cmd:_pscore}, {cmd:_weight}, {cmd:_id},
{cmd:_n1}...{cmd:_n}{it:k}, {cmd:_pdif}, {cmd:_support}, {cmd:_treated},
{cmd:_nn} and {cmd:_}{it:outcome} variables. The names {cmd:_pscore},
{cmd:_weight}, {cmd:_id}, {cmd:_support}, {cmd:_treated}, {cmd:_nn} and
{cmd:_pdif} cannot be used as inputs. Observation order is preserved.

{pstd}
{cmd:_id} and {cmd:_n}{it:k} follow {cmd:psmatch2}: {cmd:_id} is each
observation's position in the matching order (observations on support
first, then by {cmd:_treated}, by propensity score (descending with
{opt descending}), and by observation; on the Mahalanobis metric without the
score), and {cmd:_n}{it:k} is the {cmd:_id} of an observation's {it:k}th
match. They are not observation numbers: to look up the first match, sort by
{cmd:_id}, as in {cmd:sort _id} followed by
{cmd:generate double match_y = y[_n1]}.

{pstd}
When propensity scores tie (for example, with only discrete covariates),
nearest-neighbor results depend on the sort order of the data, exactly as with
{cmd:psmatch2}: among tied controls, the first in the current order is used.
To randomize the order used to break ties reproducibly, set the seed and sort
the data randomly before matching, for example {cmd:set seed 1},
{cmd:generate double u = runiform()}, {cmd:sort u}.

{pstd}
The ATT is the mean difference between the supported treated observations'
outcomes and their matched outcomes {cmd:_}{it:outcome}; with {opt ate}, the
ATU is the mean difference between the controls' matched outcomes and their
own. Observation weights are not supported.

{pstd}
Without {opt ai()}, {cmd:r(att_se)} (also returned as {cmd:r(seatt)}) uses the
default {cmd:psmatch2} ATT standard-error formula:
sqrt(s1^2/n1 + s0^2*sum(w0^2)/n1^2). Here n1 is the number of treated
observations on support, s1^2 and s0^2 are the unweighted sample outcome
variances among treated observations on support and matched controls,
and w0 is each control's matching weight. Squared weights account for control
reuse and fractional matches. This approximation treats propensity scores as
fixed and does not incorporate their estimation uncertainty. As in
{cmd:psmatch2}, it is missing with {opt spline} and with {opt llr} and the
tricube kernel, and the squared weights and the {opt ai()} terms are computed
in the default storage type ({helpb set type}).

{pstd}
The Unmatched row uses the pooled-variance standard error from an ordinary
least-squares regression of the outcome on treatment, as in {cmd:psmatch2}.


{marker examples}{...}
{title:Examples}

{pstd}Load example data and estimate treatment effect:{p_end}
{phang2}{cmd:. webuse cattaneo2, clear}{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, outcome(bweight)}{p_end}

{pstd}Use logit instead of probit for propensity score:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, outcome(bweight) logit}{p_end}

{pstd}Nearest neighbor matching with 5 neighbors:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, neighbor(5) outcome(bweight)}{p_end}

{pstd}Caliper matching with caliper = 0.05:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, caliper(0.05) outcome(bweight)}{p_end}

{pstd}ATT, ATU and ATE with Abadie-Imbens standard errors:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke, mahalanobis(mage medu) outcome(bweight) ate ai(2)}{p_end}

{pstd}Radius matching with common support:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, radius caliper(0.1) common}{p_end}

{pstd}Kernel matching with Gaussian kernel:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, kernel kerneltype(normal) bwidth(0.06)}{p_end}

{pstd}Local linear regression matching, trimming 5 percent:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, llr outcome(bweight) trim(5)}{p_end}

{pstd}Use pre-estimated propensity score:{p_end}
{phang2}{cmd:. logit mbsmoke mage medu}{p_end}
{phang2}{cmd:. predict pscore, pr}{p_end}
{phang2}{cmd:. cpsmatch mbsmoke, pscore(pscore) outcome(bweight)}{p_end}

{pstd}Matching without replacement:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, noreplacement outcome(bweight)}{p_end}

{pstd}Display timing breakdown:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, outcome(bweight) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cpsmatch} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(n_treated)}}number of treated observations{p_end}
{synopt:{cmd:r(n_controls)}}number of control observations{p_end}
{synopt:{cmd:r(n_matched)}}number of treated observations on support{p_end}
{synopt:{cmd:r(n_off_support)}}number of treated observations off support{p_end}
{synopt:{cmd:r(common_min)}}minimum control propensity score{p_end}
{synopt:{cmd:r(common_max)}}maximum control propensity score{p_end}
{synopt:{cmd:r(att)}}ATT (if outcome specified; missing when no treated observation is matched){p_end}
{synopt:{cmd:r(att_se)}}ATT standard error (see Remarks){p_end}
{synopt:{cmd:r(seatt)}}alias for {cmd:r(att_se)}, as in {cmd:psmatch2}{p_end}
{synopt:{cmd:r(att_t)}}t-statistic for ATT (if outcome specified){p_end}
{synopt:{cmd:r(atu)}, {cmd:r(ate)}}ATU and ATE (with {opt ate}){p_end}
{synopt:{cmd:r(seatu)}, {cmd:r(seate)}}their standard errors (with {opt ate}; missing without {opt ai()}){p_end}
{synopt:{cmd:r(att_}{it:outcome}{cmd:)}, {cmd:r(seatt_}{it:outcome}{cmd:)}}each outcome's ATT and standard error, as in {cmd:psmatch2} (when the name fits in 32 characters){p_end}
{synopt:{cmd:r(atu_}{it:outcome}{cmd:)}, {cmd:r(ate_}{it:outcome}{cmd:)}}each outcome's ATU and ATE (with {opt ate}){p_end}
{synopt:{cmd:r(seatu_}{it:outcome}{cmd:)}, {cmd:r(seate_}{it:outcome}{cmd:)}}their standard errors (with {opt ate}){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(method)}}matching method: nearest, radius, kernel, llr, or spline{p_end}
{synopt:{cmd:r(metric)}}pscore or mahalanobis{p_end}
{synopt:{cmd:r(treatvar)}}name of treatment variable{p_end}
{synopt:{cmd:r(outcome)}}outcome variables (if specified){p_end}
{synopt:{cmd:r(depvar)}}outcome variables, as {cmd:psmatch2} names them{p_end}
{synopt:{cmd:r(exog)}}score covariates, as in {cmd:psmatch2}{p_end}

{pstd}
{cmd:cpsmatch} creates the following variables, with {cmd:psmatch2}'s
variable labels (with the prefix {cmd:cpsmatch:}):

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Variables}{p_end}
{synopt:{cmd:_pscore}}estimated or provided propensity score (the log odds with {opt odds}){p_end}
{synopt:{cmd:_treated}}treatment indicator, value-labeled Untreated/Treated{p_end}
{synopt:{cmd:_support}}common support indicator, value-labeled Off support/On support{p_end}
{synopt:{cmd:_weight}}weight of the matched controls (with {opt ate} also of the treated): with nearest-neighbor matching, the frequency of use as a match (1/{it:k} each with {it:k} neighbors); with kernel and llr matching, the sum of the kernel weights; 1 for the treated without {opt ate}; missing for unused controls and observations off support. {cmd:char _weight[Type]} is {cmd:aweight}.{p_end}
{synopt:{cmd:_}{it:outcome}}matched outcome of each treated observation (with {opt ate}, also control){p_end}
{synopt:{cmd:_id}}position in the matching order (nearest-neighbor matching; see Remarks){p_end}
{synopt:{cmd:_n1}...{cmd:_n}{it:k}}{cmd:_id} of the {it:k} nearest neighbors{p_end}
{synopt:{cmd:_nn}}number of matches (with {opt ties} possibly more than {opt neighbor()}); zero for unmatched treated and for controls on support without {opt ate}; for radius and kernel matching on the score, the number of controls with a positive weight{p_end}
{synopt:{cmd:_pdif}}absolute score difference to the nearest neighbor (nearest-neighbor matching on the score){p_end}
{synopt:{cmd:_s_}{it:outcome}}smoothed outcome ({opt llr} with the Epanechnikov kernel, {opt spline}){p_end}
{synopt:{cmd:_self_}{it:outcome}}mean same-arm matched outcome, or the conditional variance with {opt altvariance} ({opt ai()}){p_end}


{title:References}

{phang}
Abadie, A., D. Drukker, J. L. Herr, and G. W. Imbens. 2004. Implementing
matching estimators for average treatment effects in Stata.
{it:Stata Journal} 4: 290-311.

{phang}
Abadie, A., and G. W. Imbens. 2006. Large sample properties of matching
estimators for average treatment effects. {it:Econometrica} 74: 235-267.

{phang}
Leuven, E., and B. Sianesi. 2003. PSMATCH2: Stata module to perform full
Mahalanobis and propensity score matching, common support graphing, and
covariate imbalance testing. Statistical Software Components S432001.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Online: {help psmatch2}, {help teffects psmatch}, {help logit}, {help probit}, {help ctools}
{p_end}

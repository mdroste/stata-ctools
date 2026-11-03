{smcl}
{* *! version 1.2.4 20261103}{...}
{viewerjumpto "Syntax" "cpsmatch##syntax"}{...}
{viewerjumpto "Description" "cpsmatch##description"}{...}
{viewerjumpto "Options" "cpsmatch##options"}{...}
{viewerjumpto "Remarks" "cpsmatch##remarks"}{...}
{viewerjumpto "Examples" "cpsmatch##examples"}{...}
{viewerjumpto "Stored results" "cpsmatch##results"}{...}
{title:Title}

{phang}
{bf:cpsmatch} {hline 2} C-accelerated propensity score matching for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cpsmatch}
{it:treatvar}
[{it:varlist}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Propensity Score}
{synopt:{opth ps:core(varname)}}use existing propensity score variable{p_end}
{synopt:{opt logit}}estimate propensity score using logit{p_end}
{synopt:{opt probit}}estimate propensity score using probit (default){p_end}

{syntab:Matching Method}
{synopt:{opth n:eighbor(#)}}number of nearest neighbors (default: 1){p_end}
{synopt:{opth cal:iper(#)}}caliper width for matching{p_end}
{synopt:{opt radius}}radius matching{p_end}
{synopt:{opt kernel}}kernel matching{p_end}
{synopt:{opth kerneltype(string)}}kernel function: epan, normal, biweight, uniform, tricube{p_end}
{synopt:{opth bw:idth(#)}}bandwidth for kernel matching (default: 0.06){p_end}

{syntab:Matching Options}
{synopt:{opt com:mon}}impose common support{p_end}
{synopt:{opt norep:lacement}}match without replacement{p_end}
{synopt:{opt ties}}include all tied matches{p_end}
{synopt:{opt desc:ending}}match in descending propensity score order{p_end}

{syntab:Outcome Analysis}
{synopt:{opth out:come(varname)}}outcome variable for ATT calculation{p_end}

{syntab:Performance}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cpsmatch} matches treated observations to controls using an estimated or
supplied propensity score. With {opt outcome()}, it reports the average
treatment effect on the treated (ATT) for the matched sample. Observations
with {cmd:_treated} equal to 1 are matched to those with {cmd:_treated} equal
to 0; see {help cpsmatch##remarks:Remarks} for other treatment codes and
samples without treated or control observations. Supply covariates for score
estimation or specify {opt pscore()}, not both. Covariates may include factor
variables ({help fvvarlist}).


{pstd}
Matching methods:

{p 8 12 2}1. {bf:Nearest neighbor matching:} Each treated observation is matched to
the closest control(s) based on propensity score distance.{p_end}

{p 8 12 2}2. {bf:Radius/caliper matching:} Each treated observation is matched to
all controls within a specified distance (caliper).{p_end}

{p 8 12 2}3. {bf:Kernel matching:} Controls are weighted using a kernel function
based on their propensity score distance from the treated observation.{p_end}

{pstd}
{cmd:cpsmatch} creates several variables in your dataset:

{p2colset 8 20 22 2}{...}
{p2col:{bf:_pscore}}the propensity score{p_end}
{p2col:{bf:_weight}}matching weights for each observation; missing for unused controls and observations off support{p_end}
{p2col:{bf:_id}}original observation number of the first matched control{p_end}
{p2col:{bf:_support}}common support indicator (1 = on support){p_end}
{p2col:{bf:_treated}}copy of treatment indicator{p_end}
{p2col:{bf:_nn}}number of neighbors used{p_end}


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

{dlgtab:Matching Method}

{phang}
{opth neighbor(#)} specifies the number of nearest neighbors to match to each
treated observation. The default is 1. With {opt ties}, additional matches
may be included if there are ties at the boundary.

{phang}
{opth caliper(#)} specifies the maximum propensity score distance for a
valid match. Observations without matches within the caliper are left
unmatched. The caliper is specified in absolute units of the propensity
score (e.g., {cmd:caliper(0.05)} means a maximum distance of 0.05 on the
propensity score scale). As in {cmd:psmatch2}, a match needs a distance
strictly below the caliper, and kernel matching ignores {opt caliper()}.

{phang}
{opt radius} specifies radius matching, where each treated observation is
matched to all controls within the specified caliper distance. Equal weights
are assigned to all matches. The default radius is 0.25 when {opt caliper()}
is not specified or is 0; as in {cmd:psmatch2}, a negative {opt caliper()}
matches nothing. With both {opt radius} and {opt kernel}, radius matching is
performed, as in {cmd:psmatch2}.

{phang}
{opt kernel} specifies kernel matching, where all control observations
receive a weight based on a kernel function of their propensity score
distance from the treated observation. Compact-support kernels assign zero weight outside
the bandwidth, and match nothing when {opt bwidth()} is 0 or negative; the
normal kernel has unbounded support and uses the absolute value of
{opt bwidth()}, as in {cmd:psmatch2}.

{phang}
{opth kerneltype(string)} specifies the kernel function for kernel matching:

{p 12 16 2}{bf:epan} - Epanechnikov kernel (default){p_end}
{p 12 16 2}{bf:normal} - Gaussian/normal kernel{p_end}
{p 12 16 2}{bf:biweight} - Biweight kernel{p_end}
{p 12 16 2}{bf:uniform} - Uniform kernel{p_end}
{p 12 16 2}{bf:tricube} - Tricube kernel{p_end}

{phang}
{opth bwidth(#)} specifies the bandwidth for kernel matching. The default
is 0.06. Larger bandwidths produce smoother estimates but may introduce bias.

{dlgtab:Matching Options}

{phang}
{opt common} imposes the common support restriction, excluding treated
observations with propensity scores outside the range of control propensity
scores. Controls remain on support. Treated observations without an eligible
match also receive {cmd:_support=0}.

{phang}
{opt noreplacement} specifies matching without replacement. Each control can
only be matched to one treated observation. The default is matching with
replacement. Without replacement, only nearest-neighbor matching with
{cmd:neighbor(1)} is supported. Treated observations are matched sequentially
in propensity-score order (see {opt descending}). For each one, the search
starts at the next unused control after the previous match and continues until
the score distance increases. After reaching the last control, the search
reverses direction. Treated observations remaining when no further control can
be selected receive {cmd:_support=0}. With replacement or without, the nearest
control is the first of equally distant controls in score order (the lower
score on equal distance), as in {cmd:psmatch2}.

{phang}
{opt ties} specifies that when there are tied matches at the boundary
distance, all tied observations should be included as matches.

{phang}
{opt descending} specifies that matching should proceed in descending
order of propensity score. As in {cmd:psmatch2}, it requires
{opt noreplacement}.

{dlgtab:Outcome Analysis}

{phang}
{opth outcome(varname)} specifies the outcome variable for calculating
the Average Treatment Effect on the Treated (ATT). If specified,
{cmd:cpsmatch} computes and displays the ATT estimate.

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
and, when specified, outcomes. Estimated scores also restrict matching to the
score model's estimation sample.

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
{cmd:_support}, {cmd:_treated}, and {cmd:_nn} variables. These names cannot be
used as inputs. Observation order is preserved.

{pstd}
When propensity scores tie (for example, with only discrete covariates),
nearest-neighbor results depend on the sort order of the data, exactly as with
{cmd:psmatch2}: among tied controls, the first in the current order is used.
To randomize the order used to break ties reproducibly, set the seed and sort
the data randomly before matching, for example {cmd:set seed 1},
{cmd:generate double u = runiform()}, {cmd:sort u}.


{pstd}
The ATT (Average Treatment Effect on the Treated) is calculated as the
weighted mean difference in outcomes between matched treated and control
observations. Only ATT is reported; ATE/ATU and observation-weight syntax are
not supported.

{pstd}
{cmd:r(att_se)} (also returned as {cmd:r(seatt)}) uses the default
{cmd:psmatch2} ATT standard-error formula:
sqrt(s1^2/n1 + s0^2*sum(w0^2)/n1^2). Here n1 is the number of treated
observations on support, s1^2 and s0^2 are the unweighted sample outcome
variances among treated observations on support and matched controls,
and w0 is each control's matching weight. Squared weights account for control
reuse and fractional matches. This approximation treats propensity scores as
fixed and does not incorporate their estimation uncertainty. The optional
{cmd:psmatch2} {cmd:ai()} variance estimator is not supported.

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

{pstd}Radius matching with common support:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, radius caliper(0.1) common}{p_end}

{pstd}Kernel matching with Gaussian kernel:{p_end}
{phang2}{cmd:. cpsmatch mbsmoke mage medu, kernel kerneltype(normal) bwidth(0.06)}{p_end}

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
{synopt:{cmd:r(n_matched)}}number of treated observations matched{p_end}
{synopt:{cmd:r(n_off_support)}}number of observations off common support{p_end}
{synopt:{cmd:r(common_min)}}minimum propensity score in common support{p_end}
{synopt:{cmd:r(common_max)}}maximum propensity score in common support{p_end}
{synopt:{cmd:r(att)}}estimated ATT (if outcome specified; missing when no treated observation is matched){p_end}
{synopt:{cmd:r(att_se)}}approximate ATT standard error (see Remarks){p_end}
{synopt:{cmd:r(seatt)}}alias for {cmd:r(att_se)}, as in {cmd:psmatch2}{p_end}
{synopt:{cmd:r(att_t)}}t-statistic for ATT (if outcome specified){p_end}
{synopt:{cmd:r(att_}{it:outcome}{cmd:)}}alias for {cmd:r(att)}, as in {cmd:psmatch2} (when the name fits in 32 characters){p_end}
{synopt:{cmd:r(seatt_}{it:outcome}{cmd:)}}alias for {cmd:r(att_se)}, as in {cmd:psmatch2} (when the name fits in 32 characters){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(method)}}matching method: nearest, radius, or kernel{p_end}
{synopt:{cmd:r(treatvar)}}name of treatment variable{p_end}
{synopt:{cmd:r(outcome)}}name of outcome variable (if specified){p_end}

{pstd}
{cmd:cpsmatch} creates the following variables:

{synoptset 12 tabbed}{...}
{p2col 5 12 16 2: Variables}{p_end}
{synopt:{cmd:_pscore}}estimated or provided propensity score{p_end}
{synopt:{cmd:_weight}}matching weight for each observation; missing for unused controls and observations off support, as with {cmd:psmatch2}{p_end}
{synopt:{cmd:_id}}original observation number of the first matched control (treated observations){p_end}
{synopt:{cmd:_support}}common support indicator (1 = on support, 0 = off){p_end}
{synopt:{cmd:_treated}}copy of treatment indicator{p_end}
{synopt:{cmd:_nn}}number of neighbors used for matching{p_end}


{pstd}
{cmd:_nn} records the actual number of controls used for each treated
observation. It includes ties; for radius/kernel matching it counts controls
with positive matching weight. Unmatched treated observations have zero;
controls and rows outside the matching sample have missing values.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Online: {help psmatch2}, {help teffects psmatch}, {help logit}, {help probit}, {help ctools}
{p_end}

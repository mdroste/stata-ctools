{smcl}
{* *! version 0.9.1 06Feb2026}{...}
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
{synopt:{opt v:erbose}}display detailed progress information and timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cpsmatch} performs propensity score matching to estimate causal treatment
effects in observational studies. It implements several matching methods:

{pstd}
{cmd:cpsmatch} is a high-performance replacement for {cmd:psmatch2} (v4.0.12) by
Edwin Leuven and Barbara Sianesi, available from SSC.

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
{p2col:{bf:_weight}}matching weights for each observation{p_end}
{p2col:{bf:_id}}matched observation identifier{p_end}
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
propensity score scale).

{phang}
{opt radius} specifies radius matching, where each treated observation is
matched to all controls within the specified caliper distance. Equal weights
are assigned to all matches.

{phang}
{opt kernel} specifies kernel matching, where all control observations
receive a weight based on a kernel function of their propensity score
distance from the treated observation.

{phang}
{opth kerneltype(string)} specifies the kernel function for kernel matching:

{p 12 16 2}{bf:epan} - Epanechnikov kernel (default){p_end}
{p 12 16 2}{bf:normal} or {bf:gaussian} - Gaussian/normal kernel{p_end}
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
scores (and vice versa).

{phang}
{opt noreplacement} specifies matching without replacement. Each control
can only be matched to one treated observation. The default is matching
with replacement.

{phang}
{opt ties} specifies that when there are tied matches at the boundary
distance, all tied observations should be included as matches.

{phang}
{opt descending} specifies that matching should proceed in descending
order of propensity score. This can affect results when matching without
replacement.

{dlgtab:Outcome Analysis}

{phang}
{opth outcome(varname)} specifies the outcome variable for calculating
the Average Treatment Effect on the Treated (ATT). If specified,
{cmd:cpsmatch} computes and displays the ATT estimate.

{dlgtab:Performance}

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel operations. By default, {cmd:cpsmatch} uses all available CPU cores.

{phang}
{opt verbose} displays detailed progress information and timing breakdown.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cpsmatch} uses a high-performance C plugin that parallelizes the
matching process for significantly improved performance on large datasets.

{pstd}
The propensity score is the conditional probability of receiving treatment
given observed covariates. Matching on the propensity score balances the
distribution of covariates between treated and control groups, reducing
confounding bias in treatment effect estimates.

{pstd}
Key differences from {help psmatch2:psmatch2}:

{p 8 12 2}- Parallel matching using OpenMP for multi-core speedup{p_end}
{p 8 12 2}- Optimized memory layout for cache efficiency{p_end}
{p 8 12 2}- Automatic handling of large datasets{p_end}

{pstd}
The ATT (Average Treatment Effect on the Treated) is calculated as the
weighted mean difference in outcomes between matched treated and control
observations.


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
{synopt:{cmd:r(att)}}estimated ATT (if outcome specified){p_end}
{synopt:{cmd:r(att_se)}}standard error of ATT (if outcome specified){p_end}
{synopt:{cmd:r(att_t)}}t-statistic for ATT (if outcome specified){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(method)}}matching method: nearest, radius, or kernel{p_end}
{synopt:{cmd:r(treatvar)}}name of treatment variable{p_end}
{synopt:{cmd:r(outcome)}}name of outcome variable (if specified){p_end}

{pstd}
{cmd:cpsmatch} creates the following variables:

{synoptset 12 tabbed}{...}
{p2col 5 12 16 2: Variables}{p_end}
{synopt:{cmd:_pscore}}estimated or provided propensity score{p_end}
{synopt:{cmd:_weight}}matching weight for each observation{p_end}
{synopt:{cmd:_id}}matched observation identifier (for treated: ID of matched control){p_end}
{synopt:{cmd:_support}}common support indicator (1 = on support, 0 = off){p_end}
{synopt:{cmd:_treated}}copy of treatment indicator{p_end}
{synopt:{cmd:_nn}}number of neighbors used for matching{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Online: {help psmatch2}, {help teffects psmatch}, {help logit}, {help probit}, {help ctools}
{p_end}

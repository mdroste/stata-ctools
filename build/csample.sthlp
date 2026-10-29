{smcl}
{* *! version 1.2.0 20261029}{...}
{viewerjumpto "Syntax" "csample##syntax"}{...}
{viewerjumpto "Description" "csample##description"}{...}
{viewerjumpto "Options" "csample##options"}{...}
{viewerjumpto "Remarks" "csample##remarks"}{...}
{viewerjumpto "Examples" "csample##examples"}{...}
{viewerjumpto "Stored results" "csample##results"}{...}
{title:Title}

{phang}
{bf:csample} {hline 2} C-accelerated random sampling without replacement


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:csample}
{it:#}
{ifin}
[{cmd:,} {it:options}]

{p 8 17 2}
{cmdab:csample}
{ifin}
{cmd:,} {opt count(#)}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt count(#)}}draw up to {it:#} observations{p_end}
{synopt:{opt by(varlist)}}sample within groups defined by {it:varlist}{p_end}
{synoptline}
{syntab:Advanced}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csample} draws a random sample without replacement. Within the observations
selected by {cmd:if} and {cmd:in}, it drops those not drawn. Observations outside
that selection remain in the dataset.

{pstd}
The positional {it:#} specifies a percentage from 0 to 100. The number drawn
is rounded to the nearest integer, separately within each {opt by()} group.
Use {opt count(#)} to specify a positive integer count instead. If a group has
fewer than {it:#} observations, all its observations are kept.

{pstd}
Use {help set seed} to reproduce a draw with the same data, observation order,
and options. {cmd:csample} and {cmd:sample} use different random-number generators
and need not select the same observations from the same seed.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt count(#)} specifies the number of observations to draw. This may not be
combined with the {it:#} percent syntax. When used with {opt by()}, {it:#}
observations are drawn from each group, or all observations if the group is smaller.

{phang}
{opt by(varlist)} specifies that sampling should be performed separately within
each group defined by {it:varlist}. Each group will have the specified
percentage or count of observations sampled.

{dlgtab:Advanced}

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
Sampling uses the xoshiro256** generator, seeded from Stata's random-number
state, and a partial Fisher-Yates shuffle. Numeric and string grouping variables
are allowed; missing group values form groups. The retained observations keep
their original order. Weights are not accepted.


{marker examples}{...}
{title:Examples}

{pstd}Draw a 10% random sample:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csample 10}{p_end}

{pstd}Draw exactly 20 observations:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csample, count(20)}{p_end}

{pstd}Draw 50% sample within each foreign category:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csample 50, by(foreign)}{p_end}

{pstd}Draw 5 observations from each group:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csample, count(5) by(foreign)}{p_end}

{pstd}Reproducible sampling:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. csample 25}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:csample} stores the following in global scalars (accessible via
{cmd:scalar(_csample_...)}):

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:_csample_nobs}}number of observations eligible for sampling{p_end}
{synopt:{cmd:_csample_kept}}number of eligible observations kept{p_end}
{synopt:{cmd:_csample_dropped}}number of eligible observations dropped{p_end}
{synopt:{cmd:_csample_ngroups}}number of groups{p_end}
{synopt:{cmd:_csample_threads}}number of threads used{p_end}


{marker author}{...}
{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{marker seealso}{...}
{title:Also see}

{psee}
Manual: {bf:[D] sample}

{psee}
Online: {help sample}, {help bsample}, {help cbsample}, {help ctools}
{p_end}

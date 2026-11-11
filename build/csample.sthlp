{smcl}
{* *! version 1.2.6 20261111}{...}
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
{synopt:{opt c:ount}}{it:#} is a number of observations, not a percentage{p_end}
{synopt:{opt count(#)}}draw up to {it:#} observations (same as {it:#}{cmd:, count}){p_end}
{synopt:{opt by(varlist)}}sample within groups defined by {it:varlist}{p_end}
{synopt:{opt pdup:licates(#)}}accepted for {cmd:sample} compatibility; no effect{p_end}
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
The number following {cmd:csample} specifies a percentage from 0 to 100. The
number drawn is rounded to the nearest integer, separately within each
{opt by()} group, as {cmd:sample} computes it. With {opt count}
({cmd:csample} {it:#}{cmd:, count}, as in {cmd:sample}) or {opt count(#)},
{it:#} is an integer count instead. If a group has fewer than {it:#}
observations, all its observations are kept.

{pstd}
As {cmd:sample} does, {cmd:csample} reports "({it:#} observation(s) deleted)",
does nothing (and reports nothing) for 100 percent or a count of at least
{cmd:_N}, and shows {cmd:sample}'s messages and return codes for an invalid
{it:#}. A negative count is reported and then, as in {cmd:sample}, treated as
zero.

{pstd}
Use {help set seed} to reproduce a draw with the same data, observation order,
and options. {cmd:csample} and {cmd:sample} use different random-number generators
and need not select the same observations from the same seed.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt count} specifies that {it:#} is the number of observations to draw, as
in {cmd:sample}. {opt count(#)} is the same as {it:#}{cmd:, count}; it may not
be combined with a {it:#} before the comma. When used with {opt by()}, {it:#}
observations are drawn from each group, or all observations if the group is smaller.
{cmd:count(0)}, like {cmd:sample 0, count}, drops every observation in the sample.

{phang}
{opt by(varlist)}, or the {cmd:by} {it:varlist}{cmd::} prefix as with {cmd:sample}, specifies that sampling should be performed separately within
each group defined by {it:varlist}. Each group will have the specified
percentage or count of observations sampled. {cmd:strL} group variables may
hold text of any length and binary values. As with {cmd:sample}, the data are
left sorted by the {opt by()} variables; within each group the kept
observations keep their order.

{phang}
{opt pduplicates(#)} is {cmd:sample}'s undocumented option setting how many
uniform random numbers it draws per observation to avoid ties. {cmd:csample}
draws without ties, so the option is accepted and has no effect.

{dlgtab:Advanced}

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256. Zero uses the default.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
Sampling uses the xoshiro256** generator, seeded from Stata's random-number
state, and a partial Fisher-Yates shuffle. Numeric and string grouping variables
are allowed; missing group values form groups. Without {opt by()}, the retained
observations keep their original order ({cmd:sample} leaves them in random
order). Weights are not accepted.

{pstd}
Under {help version} 6 or earlier, {cmd:sample} uses its old algorithm, so that
old scripts reproduce their samples. {cmd:csample} then calls {cmd:sample}, and
the draw, messages and return codes are {cmd:sample}'s; {opt threads()} and
{opt verbose} are ignored.


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

{pstd}
An empty eligible sample publishes zero counts and timings, with zero threads.

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

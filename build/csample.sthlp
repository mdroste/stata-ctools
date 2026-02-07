{smcl}
{* *! version 0.9.1 06Feb2026}{...}
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
{synopt:{opt count(#)}}draw exactly {it:#} observations{p_end}
{synopt:{opt by(varlist)}}sample within groups defined by {it:varlist}{p_end}
{synoptline}
{syntab:Advanced}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csample} is a high-performance drop-in replacement for Stata's {help sample:sample}
command. It draws a random sample without replacement from the dataset in memory,
dropping observations not selected.

{pstd}
{cmd:csample} {it:#} draws a {it:#}% pseudorandom sample of the data, where
0 <= {it:#} <= 100. Observations not in the sample are dropped.

{pstd}
{cmd:csample, count(}{it:#}{cmd:)} draws a sample of {it:#} observations.
When combined with {opt by()}, {it:#} observations are drawn from each group.

{pstd}
{cmd:csample} respects Stata's {help set seed:set seed} for reproducibility.
Setting the same seed before calling {cmd:csample} will produce identical samples.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt count(#)} specifies the number of observations to draw. This may not be
combined with the {it:#} percent syntax. When used with {opt by()}, {it:#}
observations are drawn from each group.

{phang}
{opt by(varlist)} specifies that sampling should be performed separately within
each group defined by {it:varlist}. Each group will have the specified
percentage or count of observations sampled.

{dlgtab:Advanced}

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:csample} uses all available CPU cores.

{phang}
{opt verbose} displays detailed progress information and timing breakdown.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:csample} uses a high-quality xoshiro256** pseudorandom number generator
and Fisher-Yates partial shuffle algorithm for statistically reliable sampling.

{pstd}
When using {opt by()}, groups are processed in parallel for improved performance.

{pstd}
Unlike {help sample:sample}, {cmd:csample} does not support frequency weights.


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
{synopt:{cmd:_csample_kept}}number of observations kept{p_end}
{synopt:{cmd:_csample_dropped}}number of observations dropped{p_end}
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

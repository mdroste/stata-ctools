{smcl}
{* *! version 1.0.1 07Feb2026}{...}
{viewerjumpto "Syntax" "cbsample##syntax"}{...}
{viewerjumpto "Description" "cbsample##description"}{...}
{viewerjumpto "Options" "cbsample##options"}{...}
{viewerjumpto "Remarks" "cbsample##remarks"}{...}
{viewerjumpto "Examples" "cbsample##examples"}{...}
{viewerjumpto "Stored results" "cbsample##results"}{...}
{title:Title}

{phang}
{bf:cbsample} {hline 2} C-accelerated bootstrap sampling with replacement


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cbsample}
{ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt n(#)}}size of bootstrap sample; default is {cmd:n(_N)}{p_end}
{synopt:{opt cl:uster(varlist)}}variables identifying resampling clusters{p_end}
{synopt:{opt str:ata(varlist)}}variables identifying strata{p_end}
{synopt:{opt weight(newvar)}}store bootstrap weights instead of resampling{p_end}
{synoptline}
{syntab:Advanced}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cbsample} is a high-performance replacement for Stata's {help bsample:bsample}
command. It draws a bootstrap sample (with replacement) from the dataset in memory.

{pstd}
By default, {cmd:cbsample} replaces the data in memory with a bootstrap sample
of the same size as the original. Observations not selected are dropped, while
selected observations may appear multiple times (implemented by keeping one copy
and dropping unselected observations by default).

{pstd}
With {opt weight()}, {cmd:cbsample} keeps all observations and stores
frequency weights indicating how many times each observation was selected.
This is useful for bootstrap estimation where weights are preferred over
physical replication of observations.

{pstd}
{cmd:cbsample} respects Stata's {help set seed:set seed} for reproducibility.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt n(#)} specifies the size of the bootstrap sample. The default is {cmd:n(_N)},
meaning a sample of the same size as the original dataset.

{phang}
{opt cluster(varlist)} specifies the variables that identify resampling clusters.
When specified, entire clusters are resampled rather than individual observations.
All observations within a selected cluster receive the same weight.

{phang}
{opt strata(varlist)} specifies variables identifying strata. Resampling is
performed independently within each stratum, maintaining the relative size of
each stratum in the bootstrap sample.

{phang}
{opt weight(newvar)} creates a new variable {it:newvar} containing bootstrap
frequency weights. When this option is specified, all observations are kept and
the weight variable indicates how many times each observation was selected
(0 if not selected). This is more memory-efficient than physical replication
for bootstrap estimation.

{dlgtab:Advanced}

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:cbsample} uses all available CPU cores.

{phang}
{opt verbose} displays detailed progress information and timing breakdown.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cbsample} uses a high-quality xoshiro256** pseudorandom number generator
for statistically reliable bootstrap sampling.

{pstd}
When neither {opt cluster()} nor {opt strata()} is specified, simple random
sampling with replacement is performed.

{pstd}
With {opt cluster()}, the cluster bootstrap resamples clusters rather than
observations. This is appropriate when observations within clusters are not
independent.

{pstd}
With {opt strata()}, stratified bootstrap sampling maintains the stratum
structure of the original data.


{marker examples}{...}
{title:Examples}

{pstd}Simple bootstrap sample:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample}{p_end}

{pstd}Bootstrap sample of specific size:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample, n(50)}{p_end}

{pstd}Store bootstrap weights instead of resampling:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample, weight(bsweight)}{p_end}
{phang2}{cmd:. regress price mpg [fw=bsweight]}{p_end}

{pstd}Stratified bootstrap:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample, strata(foreign) weight(bsweight)}{p_end}

{pstd}Cluster bootstrap:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. cbsample, cluster(idcode) weight(bsweight)}{p_end}

{pstd}Reproducible bootstrap:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. cbsample}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cbsample} stores the following in global scalars (accessible via
{cmd:scalar(_cbsample_...)}):

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:_cbsample_n_selected}}number of observations in bootstrap sample{p_end}
{synopt:{cmd:_cbsample_nclusters}}number of clusters (if {opt cluster()} specified){p_end}
{synopt:{cmd:_cbsample_nstrata}}number of strata (if {opt strata()} specified){p_end}


{marker author}{...}
{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{marker seealso}{...}
{title:Also see}

{psee}
Manual: {bf:[R] bsample}

{psee}
Online: {help bsample}, {help bootstrap}, {help csample}, {help ctools}
{p_end}

{smcl}
{* *! version 1.1.0 20261016}{...}
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
[{it:#}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{it:#}}positional number of draws within each stratum; default is the original number of sampling units{p_end}
{synopt:{opt cl:uster(varlist)}}variables identifying resampling clusters{p_end}
{synopt:{opt str:ata(varlist)}}variables identifying strata{p_end}
{synopt:{opt weight(varname)}}store bootstrap weights instead of resampling{p_end}
{synoptline}
{syntab:Advanced}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cbsample} draws a sample with replacement. A sampling unit is an observation,
or a cluster when {opt cluster()} is specified. By default, the number of draws
in each stratum equals its original number of sampling units.

{pstd}
The sample replaces the data in memory. Selected observations are replicated
according to their draw counts; unselected observations, including those outside
{cmd:if}/{cmd:in} or with missing strata, are dropped. With unequal cluster
sizes, the resulting number of observations can differ from the original.

{pstd}
With {opt weight()}, all observations remain in memory and an existing numeric
variable receives each observation's draw count. Observations outside the
sampling sample receive zero.

{pstd}
Use {help set seed} to reproduce a draw with the same data, observation order,
and options. {cmd:cbsample} and {cmd:bsample} need not produce the same sample
from the same seed.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
The positional {it:#} specifies draws within each stratum. For example,
{cmd:cbsample 2, strata(s)} draws two observations from every stratum;
adding {cmd:cluster(c)} draws two clusters from every stratum. Omitting the
count uses each stratum's original number of sampling units. An explicit count
must be a positive integer no greater than the available units in any stratum.

{phang}
{opt cluster(varlist)} specifies the variables that identify resampling clusters.
When specified, entire clusters are resampled rather than individual observations.
All observations within a selected cluster receive the same weight.

{phang}
{opt strata(varlist)} specifies variables identifying strata. Resampling is
performed independently within each stratum. Missing strata are excluded.
Cluster labels reused across strata identify separate sampling units.

{phang}
{opt weight(varname)} overwrites an existing numeric variable with bootstrap
frequency weights. Create it first with {cmd:generate double bsweight = 0}. When this option is specified, all observations are kept and
the weight variable indicates how many times each observation was selected
(0 if not selected). This is more memory-efficient than physical replication
for bootstrap estimation.

{dlgtab:Advanced}

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
The command uses the xoshiro256** generator seeded from Stata's random-number
state. With {opt cluster()}, all eligible observations in a cluster receive
the same draw count. With {opt strata()}, draws are independent across strata.


{marker examples}{...}
{title:Examples}

{pstd}Simple bootstrap sample:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample}{p_end}

{pstd}Bootstrap sample of specific size:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbsample 50}{p_end}

{pstd}Store bootstrap weights instead of resampling:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. generate double bsweight = 0}{p_end}
{phang2}{cmd:. cbsample, weight(bsweight)}{p_end}
{phang2}{cmd:. regress price mpg [fw=bsweight]}{p_end}

{pstd}Stratified bootstrap:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. generate double bsweight = 0}{p_end}
{phang2}{cmd:. cbsample, strata(foreign) weight(bsweight)}{p_end}

{pstd}Cluster bootstrap:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. generate double bsweight = 0}{p_end}
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
{synopt:{cmd:_cbsample_n_selected}}number of distinct observations selected at least once{p_end}
{synopt:{cmd:_cbsample_total_weight}}sum of draw counts across observations; size of the expanded sample{p_end}
{synopt:{cmd:_cbsample_nobs}}number of eligible observations{p_end}
{synopt:{cmd:_cbsample_threads}}number of threads used{p_end}
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

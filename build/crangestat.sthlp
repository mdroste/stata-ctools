{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "crangestat##syntax"}{...}
{viewerjumpto "Description" "crangestat##description"}{...}
{viewerjumpto "Options" "crangestat##options"}{...}
{viewerjumpto "Statistics" "crangestat##statistics"}{...}
{viewerjumpto "Remarks" "crangestat##remarks"}{...}
{viewerjumpto "Examples" "crangestat##examples"}{...}
{viewerjumpto "Stored results" "crangestat##results"}{...}
{title:Title}

{phang}
{bf:crangestat} {hline 2} C-accelerated range statistics for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:crangestat}
{it:stat_specification} [{it:stat_specification} ...]
{cmd:,}
{opt int:erval(keyvar low high)}
[{it:options}]

{pstd}
where {it:stat_specification} is

{p 8 17 2}
{cmd:(}{it:statname}{cmd:)} [{it:newvar}{cmd:=}]{varname}

{pstd}
Without an explicit {it:newvar}, the result is named {it:varname_statname},
matching {cmd:rangestat}; for example, {cmd:(mean) price} creates {cmd:price_mean}.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt int:erval(keyvar low high)}}key variable and interval bounds{p_end}

{syntab:Optional}
{synopt:{opt by(varlist)}}compute statistics within groups{p_end}
{synopt:{opt excludeself}}exclude current observation from calculations{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:crangestat} calculates statistics for each observation using all observations
where a numeric key variable is within the low and high bounds defined for the
current observation. The result variables are added to the current dataset; observation order is
preserved.

{pstd}
For each observation, {cmd:crangestat} finds all observations where:

{p 8 12 2}
{it:keyvar}[current] + {it:low} <= {it:keyvar}[other] <= {it:keyvar}[current] + {it:high}

{pstd}
and computes the specified statistics using the values of the source variables
from those observations.

{pstd}
This is particularly useful for:

{p 8 12 2}- Rolling window statistics over time series{p_end}
{p 8 12 2}- Computing statistics using observations with similar values{p_end}
{p 8 12 2}- Calculating leave-one-out statistics with {opt excludeself}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt interval(keyvar low high)} specifies the key variable and the interval
bounds. The key variable must be numeric. Bounds must be numeric constants
or {cmd:.}; variable bounds are not accepted. For each observation, other
observations are included in the calculation if their key value falls within
[current_key + low, current_key + high].

{p 8 8 2}
Use {cmd:.} (missing) for {it:low} to indicate negative infinity, or for
{it:high} to indicate positive infinity. For example, {cmd:interval(year . 0)}
includes all observations from the earliest date up to the current observation.

{dlgtab:Optional}

{phang}
{opt by(varlist)} specifies that statistics should be computed separately
within groups defined by the variables in {it:varlist}. Observations in
different groups are never included in each other's calculations.

{phang}
{opt excludeself} specifies that the current observation should be excluded
from its own calculation. This is useful for computing leave-one-out
statistics like the mean of all {it:other} observations in the range.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.

{dlgtab:Data recovery}

{pstd}
By default, the command removes new result variables after an error or
interruption.

{phang}
{opt snapshot} uses {help preserve} to copy the whole dataset for recovery.
Stata may write this copy to disk when it exceeds {cmd:c(max_preservemem)}.

{phang}
{opt fast} skips recovery. An error or interruption can leave partially filled
result variables. It cannot be combined with {opt snapshot}.


{marker statistics}{...}
{title:Statistics}

{pstd}
The following statistics are available:

{synoptset 15 tabbed}{...}
{p2col :{it:statname}}Description{p_end}
{p2line}
{p2col :{opt count}}number of non-missing observations in range{p_end}
{p2col :{opt mean}}arithmetic mean{p_end}
{p2col :{opt sum}}sum{p_end}
{p2col :{opt min}}minimum{p_end}
{p2col :{opt max}}maximum{p_end}
{p2col :{opt sd}}standard deviation (N-1 denominator){p_end}
{p2col :{opt variance}}variance (N-1 denominator){p_end}
{p2col :{opt median}}median (50th percentile){p_end}
{p2col :{opt iqr}}interquartile range (p75 - p25){p_end}
{p2col :{opt first}}first value, including missing (by key variable order){p_end}
{p2col :{opt last}}last value, including missing (by key variable order){p_end}
{p2col :{opt firstnm}}first non-missing value{p_end}
{p2col :{opt lastnm}}last non-missing value{p_end}
{p2col :{opt p1}}1st percentile{p_end}
{p2col :{opt p5}}5th percentile{p_end}
{p2col :{opt p10}}10th percentile{p_end}
{p2col :{opt p25}}25th percentile{p_end}
{p2col :{opt p75}}75th percentile{p_end}
{p2col :{opt p90}}90th percentile{p_end}
{p2col :{opt p95}}95th percentile{p_end}
{p2col :{opt p99}}99th percentile{p_end}
{p2col :{opt skewness}}skewness{p_end}
{p2col :{opt kurtosis}}kurtosis (fourth standardized moment; normal distribution = 3){p_end}
{p2line}


{marker remarks}{...}
{title:Remarks}

{pstd}
Each statistic specification accepts one numeric source variable. Repeat the
parenthesized statistic for additional variables. Output names must be distinct
new variable names. The command does not accept {cmd:if}, {cmd:in}, or weights.

{pstd}
Missing interval keys are excluded from all windows and receive missing
results. Missing source values are excluded except for {opt first} and
{opt last}. These two statistics use key order, with ties in original
observation order. {opt excludeself} removes only the current observation,
not other observations with the same key.

{pstd}
A window containing no observations produces missing results. In a nonempty
window with no nonmissing source values after exclusions, {opt count} and
{opt sum} return zero. Standard
deviation and variance require two nonmissing values. Skewness and kurtosis
use population moments and are missing when the variance is zero.


{marker examples}{...}
{title:Examples}

{pstd}Rolling mean over a 5-period window:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. gen obs = _n}{p_end}
{phang2}{cmd:. crangestat (mean) price_ma=price, interval(obs -2 2)}{p_end}

{pstd}Mean of observations with similar values (within +/- 1):{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. crangestat (mean) similar_wage=ln_wage, interval(age -1 1)}{p_end}

{pstd}Leave-one-out mean (exclude current observation):{p_end}
{phang2}{cmd:. crangestat (mean) loo_wage=ln_wage, interval(age -1 1) excludeself}{p_end}

{pstd}Multiple statistics at once:{p_end}
{phang2}{cmd:. crangestat (mean) mean_wage=ln_wage (sd) sd_wage=ln_wage (count) n_wage=ln_wage, interval(year -2 0)}{p_end}

{pstd}Rolling statistics by group:{p_end}
{phang2}{cmd:. webuse grunfeld, clear}{p_end}
{phang2}{cmd:. crangestat (mean) invest_ma=invest (sum) invest_sum=invest, interval(year -2 0) by(company)}{p_end}

{pstd}Cumulative statistics (from beginning to current):{p_end}
{phang2}{cmd:. crangestat (mean) cumul_mean=invest (max) cumul_max=invest, interval(year . 0) by(company)}{p_end}

{pstd}Forward-looking statistics:{p_end}
{phang2}{cmd:. crangestat (mean) future_mean=invest, interval(year 0 2) by(company)}{p_end}

{pstd}With verbose timing output:{p_end}
{phang2}{cmd:. crangestat (mean) invest_recent=invest, interval(year -5 0) by(company) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
When the {opt verbose} option is specified, {cmd:crangestat} stores the
following in global scalars:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:_crangestat_time_load}}time to load data from Stata to C{p_end}
{synopt:{cmd:_crangestat_time_sort}}time to sort data by key variable{p_end}
{synopt:{cmd:_crangestat_time_groups}}time to detect groups{p_end}
{synopt:{cmd:_crangestat_time_compute}}time to compute statistics{p_end}
{synopt:{cmd:_crangestat_time_store}}time to store results back to Stata{p_end}
{synopt:{cmd:_crangestat_time_total}}total plugin execution time{p_end}
{synopt:{cmd:_crangestat_nobs}}number of observations processed{p_end}
{synopt:{cmd:_crangestat_nstats}}number of statistics computed{p_end}
{synopt:{cmd:_crangestat_ngroups}}number of groups{p_end}
{synopt:{cmd:_crangestat_threads}}number of threads used{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command is inspired by {cmd:rangestat} by Robert Picard, Nicholas J. Cox,
and Roberto Ferrer. The original {cmd:rangestat} package is available from SSC:
{browse "https://ideas.repec.org/c/boc/bocode/s458161.html"}.


{title:Also see}

{psee}
Online: {browse "https://ideas.repec.org/c/boc/bocode/s458161.html":rangestat} (if installed),
{help summarize}, {help egen}, {help tsegen}, {help ctools}
{p_end}

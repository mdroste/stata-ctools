{smcl}
{* *! version 1.0.1 07Feb2026}{...}
{viewerjumpto "Syntax" "cwinsor##syntax"}{...}
{viewerjumpto "Description" "cwinsor##description"}{...}
{viewerjumpto "Options" "cwinsor##options"}{...}
{viewerjumpto "Remarks" "cwinsor##remarks"}{...}
{viewerjumpto "Examples" "cwinsor##examples"}{...}
{viewerjumpto "Stored results" "cwinsor##results"}{...}
{title:Title}

{phang}
{bf:cwinsor} {hline 2} C-accelerated parallel winsorization for Stata variables


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cwinsor}
{varlist}
{ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Percentile options}
{synopt:{opt cuts(# #)}}lower and upper percentile cutoffs{p_end}
{synopt:{opt p(#)}}lower percentile cutoff (default 1){p_end}
{synopt:{opt q(#)}}upper percentile cutoff (default 99){p_end}

{syntab:Treatment options}
{synopt:{opt trim}}trim (set to missing) instead of winsorize{p_end}

{syntab:By-group options}
{synopt:{opt by(varlist)}}winsorize within groups defined by {it:varlist}{p_end}

{syntab:Output options}
{synopt:{opt suf:fix(string)}}generate new variables with suffix{p_end}
{synopt:{opt pre:fix(string)}}generate new variables with prefix{p_end}
{synopt:{opt replace}}replace existing variables (default){p_end}

{syntab:Performance options}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cwinsor} is a high-performance C-accelerated replacement for {cmd:winsor2} (v1.1) by Yujun Lian, available from SSC.
It winsorizes or trims numeric variables at specified percentiles using parallel
processing and an efficient O(n) percentile algorithm.

{pstd}
{bf:Winsorization} limits extreme values by replacing observations below the lower
percentile with the lower percentile value, and observations above the upper
percentile with the upper percentile value. This is useful for reducing the
influence of outliers in statistical analysis.

{pstd}
{bf:Trimming} (with the {opt trim} option) instead sets extreme values to missing,
effectively removing them from subsequent analyses.


{marker options}{...}
{title:Options}

{dlgtab:Percentile options}

{phang}
{opt cuts(# #)} specifies the lower and upper percentile cutoffs as two numbers.
For example, {cmd:cuts(5 95)} winsorizes at the 5th and 95th percentiles. This
option takes precedence over {opt p()} and {opt q()}.

{phang}
{opt p(#)} specifies the lower percentile cutoff. The default is 1, meaning values
below the 1st percentile are winsorized.

{phang}
{opt q(#)} specifies the upper percentile cutoff. The default is 99, meaning values
above the 99th percentile are winsorized.

{dlgtab:Treatment options}

{phang}
{opt trim} specifies that extreme values should be trimmed (set to missing) rather
than winsorized (clamped to the percentile bounds). This is useful when you want
to completely remove outliers from analysis rather than reduce their influence.

{dlgtab:By-group options}

{phang}
{opt by(varlist)} specifies that winsorization should be performed separately
within groups defined by the variables in {it:varlist}. Percentile cutoffs are
computed independently for each group. This is essential when your data contains
multiple categories or panels where the distribution differs across groups.

{dlgtab:Output options}

{phang}
{opt suffix(string)} generates new variables with the specified suffix instead of
replacing the original variables. For example, {cmd:suffix(_w)} creates new
variables {it:varname}_w containing the winsorized values.

{phang}
{opt prefix(string)} generates new variables with the specified prefix instead of
replacing the original variables. For example, {cmd:prefix(w_)} creates new
variables w_{it:varname} containing the winsorized values.

{phang}
{opt replace} replaces the original variables with winsorized values. This is the
default behavior when neither {opt suffix()} nor {opt prefix()} is specified.
Cannot be combined with {opt suffix()} or {opt prefix()}.

{dlgtab:Performance options}

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:cwinsor} uses all available CPU cores.

{phang}
{opt verbose} displays detailed progress information and timing breakdown.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cwinsor} achieves high performance through several optimizations:

{p 8 12 2}1. {bf:Parallel processing}: Multiple variables are winsorized
simultaneously using OpenMP parallelization.{p_end}

{p 8 12 2}2. {bf:O(n) percentile computation}: Uses the quickselect algorithm
to compute percentiles in linear time on average, avoiding the O(n log n)
cost of full sorting.{p_end}

{p 8 12 2}3. {bf:Efficient memory access}: Cache-friendly data layout and
SIMD-optimized loops for maximum throughput.{p_end}

{pstd}
For large datasets (millions of observations), {cmd:cwinsor} can be 5-20x
faster than {cmd:winsor2}, depending on the number of variables, groups,
and available CPU cores.

{pstd}
{bf:Missing values} are handled correctly: they are excluded from percentile
calculations but preserved in the output (unless the observation is trimmed
due to being an outlier).

{pstd}
{bf:Minimum observations}: Winsorization requires at least 3 non-missing
observations in each group. Groups with fewer observations are left unchanged.


{marker examples}{...}
{title:Examples}

{pstd}Basic winsorization at 1st and 99th percentiles (default):{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cwinsor price}{p_end}

{pstd}Winsorize at 5th and 95th percentiles:{p_end}
{phang2}{cmd:. cwinsor price mpg, cuts(5 95)}{p_end}

{pstd}Equivalent using p() and q():{p_end}
{phang2}{cmd:. cwinsor price mpg, p(5) q(95)}{p_end}

{pstd}Trim instead of winsorize (set outliers to missing):{p_end}
{phang2}{cmd:. cwinsor price, trim}{p_end}

{pstd}Winsorize by groups:{p_end}
{phang2}{cmd:. cwinsor price, by(foreign)}{p_end}

{pstd}Generate new variables instead of replacing:{p_end}
{phang2}{cmd:. cwinsor price mpg, suffix(_w)}{p_end}
{phang2}{cmd:. list price price_w mpg mpg_w in 1/10}{p_end}

{pstd}Generate with prefix:{p_end}
{phang2}{cmd:. cwinsor price, prefix(w_)}{p_end}

{pstd}Winsorize with verbose timing output:{p_end}
{phang2}{cmd:. cwinsor price mpg weight, verbose}{p_end}

{pstd}Winsorize by groups at 2.5th and 97.5th percentiles:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. cwinsor ln_wage, by(race) cuts(2.5 97.5)}{p_end}

{pstd}Limit parallelism to 4 threads:{p_end}
{phang2}{cmd:. cwinsor price, threads(4)}{p_end}

{pstd}Winsorize multiple variables by group with new variable creation:{p_end}
{phang2}{cmd:. cwinsor price mpg weight, by(foreign) suffix(_w) cuts(1 99)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cwinsor} does not store results in {cmd:r()} or {cmd:e()}.

{pstd}
When the {opt verbose} option is specified, timing information is displayed
but not stored.


{title:Comparison with winsor2}

{pstd}
{cmd:cwinsor} provides the same functionality as {cmd:winsor2} but with
significantly better performance:

{p2colset 5 25 27 2}{...}
{p2col :{it:Feature}}{it:cwinsor}{space 8}{it:winsor2}{p_end}
{p2line}
{p2col :Parallel processing}Yes{space 10}No{p_end}
{p2col :Percentile algorithm}O(n) quickselect{space 2}O(n log n) sort{p_end}
{p2col :By-group support}Yes{space 10}Yes{p_end}
{p2col :Trim option}Yes{space 10}Yes{p_end}
{p2col :Generate new vars}Yes{space 10}Yes{p_end}
{p2col :Speed (large data)}5-20x faster{space 4}Baseline{p_end}


{title:Technical details}

{pstd}
The quickselect algorithm used for percentile computation has O(n) average-case
time complexity, compared to O(n log n) for sort-based approaches. This makes
a significant difference when processing many variables or groups, as each
percentile computation is much faster.

{pstd}
The implementation uses median-of-three pivot selection to avoid worst-case
O(n^2) behavior, and switches to insertion sort for small subarrays for
optimal performance.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command provides similar functionality to {cmd:winsor2} by Yujun Lian.
The original {cmd:winsor2} package is available from SSC:
{browse "https://ideas.repec.org/c/boc/bocode/s457765.html"}. We thank Yujun Lian
for his useful contribution to the Stata community.


{title:Also see}

{psee}
Online: {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2} (if installed),
{help summarize}, {help centile}, {help pctile}, {help ctools}
{p_end}

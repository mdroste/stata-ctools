{smcl}
{* *! version 1.1.0 20261016}{...}
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

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

{syntab:Performance options}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cwinsor} winsorizes numeric variables at specified percentiles. Values below
the lower cutoff are replaced with that cutoff; values above the upper cutoff
are replaced with the upper cutoff. With {opt trim}, those values become missing.
The default cutoffs are the 1st and 99th percentiles.

{pstd}
Without {opt prefix()} or {opt suffix()}, the command changes the source
variables. New variables inherit the source variable label, or its name when
unlabeled. New variables are missing outside {cmd:if}/{cmd:in}; in-place changes
leave observations outside that sample unchanged.


{marker options}{...}
{title:Options}

{dlgtab:Percentile options}

{phang}
{opt cuts(# #)} specifies the lower and upper percentile cutoffs as two numbers.
For example, {cmd:cuts(5 95)} winsorizes at the 5th and 95th percentiles. This
option takes precedence over {opt p()} and {opt q()}. Cutoffs must lie between
0 and 100, with the lower strictly less than the upper. Use 0 or 100 to leave
that tail unchanged.

{phang}
{opt p(#)} specifies the lower percentile cutoff. The default is 1, meaning values
below the 1st percentile are winsorized.

{phang}
{opt q(#)} specifies the upper percentile cutoff. The default is 99, meaning values
above the 99th percentile are winsorized.

{dlgtab:Treatment options}

{phang}
{opt trim} specifies that extreme values should be trimmed (set to missing) rather
than winsorized (clamped to the percentile bounds).

{dlgtab:By-group options}

{phang}
{opt by(varlist)} specifies that winsorization should be performed separately
within groups defined by the variables in {it:varlist}. Percentile cutoffs are
computed independently for each group.

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

{dlgtab:Data recovery}

{pstd}
By default, the command restores the data after an error or interruption.
In-place operations retain a copy of the values being changed until completion.

{phang}
{opt snapshot} uses {help preserve} to copy the whole dataset for recovery.
Stata may write this copy to disk when it exceeds {cmd:c(max_preservemem)}.

{phang}
{opt fast} skips recovery to reduce memory use. An error or interruption can
leave changed values or partially filled new variables. It cannot be combined
with {opt snapshot}.


{dlgtab:Performance options}

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
Percentiles are computed from nonmissing observations, separately for each
variable and, with {opt by()}, each group. Groups with one or two nonmissing
observations are included. Missing numeric or string group values form groups.

{pstd}
The percentile rule follows {cmd:winsor2}. For sorted nonmissing values
x[1], ..., x[N] and p between 0 and 100, let h=N*p/100. If h is an integer,
the cutoff is (x[h]+x[h+1])/2; otherwise it is x[ceil(h)]. The endpoints
are the minimum and maximum. Ungrouped cutoffs pass through Stata numeric
local macros, which can round values at a cutoff.

{pstd}
New variables from grouped winsorization use Stata's default numeric storage
type. Other outputs retain the source storage type. Missing source values do not
enter percentile calculations; generated outputs use system missing for them.


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
{phang2}{cmd:. sysuse auto, clear}{p_end}
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


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command provides similar functionality to {cmd:winsor2} by Yujun Lian.
The original {cmd:winsor2} package is available from SSC:
{browse "https://ideas.repec.org/c/boc/bocode/s457765.html"}.


{title:Also see}

{psee}
Online: {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2} (if installed),
{help summarize}, {help centile}, {help pctile}, {help ctools}
{p_end}

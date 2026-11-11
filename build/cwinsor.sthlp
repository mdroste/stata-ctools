{smcl}
{* *! version 1.2.6 20261111}{...}
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
{weight}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Percentile options}
{synopt:{opt c:uts(# #)}}lower and upper percentile cutoffs{p_end}
{synopt:{opt p(#)}}lower percentile cutoff (default 1){p_end}
{synopt:{opt q(#)}}upper percentile cutoff (default 99){p_end}

{syntab:Treatment options}
{synopt:{opt t:rim}}trim (set to missing) instead of winsorize{p_end}

{syntab:By-group options}
{synopt:{opt by(varlist)}}winsorize within groups defined by {it:varlist}{p_end}
{synopt:{opt nomiss:ing}}with {opt by()}, leave out groups with a missing by-value{p_end}

{syntab:Output options}
{synopt:{opt s:uffix(string)}}generate new variables with suffix{p_end}
{synopt:{opt p:refix(string)}}generate new variables with prefix{p_end}
{synopt:{opt gen:erate(namelist)}}generate the named variables, one per source{p_end}
{synopt:{opt replace}}replace existing variables (default); with new variables, overwrite existing ones{p_end}
{synopt:{opt l:abel}}add the cuts to the variable labels, as {cmd:winsor2} does{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

{syntab:Performance options}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}

{syntab:gtools options (accepted for {cmd:gstats winsor} scripts)}
{synopt:{opt bench:mark}, {opt bench:marklevel(#)}}display timing breakdown, as {opt verbose}{p_end}
{synopt:{opt compress}, {opt forcestrl}}no effect{p_end}
{synopt:{opt hash:method(str)}, {opt oncollision(str)}}checked, no effect{p_end}
{synoptline}
{p 4 6 2}
{opt aweight}s, {opt fweight}s, {opt iweight}s, and {opt pweight}s are allowed,
as in {cmd:gstats winsor}; see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cwinsor} winsorizes numeric variables at specified percentiles. Values below
the lower cutoff are replaced with that cutoff; values above the upper cutoff
are replaced with the upper cutoff. With {opt trim}, those values become missing.
The default cutoffs are the 1st and 99th percentiles.

{pstd}
Without {opt prefix()}, {opt suffix()} or {opt generate()}, the command changes
the source variables. New variables inherit the source variable label, or its name when
unlabeled. New variables are missing outside {cmd:if}/{cmd:in}; in-place changes
leave observations outside that sample unchanged.


{marker options}{...}
{title:Options}

{dlgtab:Percentile options}

{phang}
{opt cuts(# #)} specifies the lower and upper percentile cutoffs as two numbers.
For example, {cmd:cuts(5 95)} winsorizes at the 5th and 95th percentiles. This
option takes precedence over {opt p()} and {opt q()}. Cutoffs must lie between
0 and 100, the lower first. Use 0 or 100 to leave that tail unchanged. A
reversed pair such as {cmd:cuts(99 1)} exits with r(198), as in
{cmd:gstats winsor} ({cmd:winsor2} meant to swap it, but stops with r(111)).

{phang}
{opt p(#)} specifies the lower percentile cutoff. The default is 1, meaning values
below the 1st percentile are winsorized. {opt p()} is also the abbreviation of
{opt prefix()}, as in {cmd:gstats winsor}: a number is the percentile, a name
the prefix.

{phang}
{opt q(#)} specifies the upper percentile cutoff. The default is 99, meaning values
above the 99th percentile are winsorized. As in {cmd:winsor2}, it may equal, but not
be less than, the lower cutoff.

{dlgtab:Treatment options}

{phang}
{opt trim} sets values outside the percentile cutoffs to missing instead of
replacing them with the cutoffs.

{dlgtab:By-group options}

{phang}
{opt by(varlist)} specifies that winsorization should be performed separately
within groups defined by the variables in {it:varlist}. Percentile cutoffs are
computed independently for each group.

{phang}
{opt nomissing}, with {opt by()}, leaves out the observations with a missing
value of any by-variable, as {cmd:gstats winsor} does: they are neither used
nor changed (new variables are missing there). Without it, missing by-values
form groups of their own.

{dlgtab:Output options}

{phang}
{opt suffix(string)} generates new variables with the specified suffix instead of
replacing the original variables. For example, {cmd:suffix(_w)} creates new
variables {it:varname}_w containing the winsorized values.

{phang}
{opt prefix(string)} generates new variables with the specified prefix instead of
replacing the original variables. For example, {cmd:prefix(w_)} creates new
variables w_{it:varname} containing the winsorized values. With {opt suffix()}
too, as in {cmd:gstats winsor}, the new variables are
{it:prefix}{it:varname}{it:suffix}.

{pstd}
The display format of a new variable follows the reference command that has
the option used: with {opt prefix()}, {opt generate()} or weights, which only
{cmd:gstats winsor} has, it is the source's format, as there; otherwise it is
the default format of its storage type, as {cmd:winsor2}'s {cmd:generate}
gives.

{phang}
{opt generate(namelist)} generates the named variables, the {it:i}th holding
the winsorized {it:i}th variable of {it:varlist}, as in {cmd:gstats winsor}.
It takes as many names as {it:varlist} has variables and cannot be combined with
{opt suffix()} or {opt prefix()}.

{phang}
{opt replace} replaces the original variables with winsorized values. This is the
default behavior when none of {opt suffix()}, {opt prefix()} and {opt generate()}
is specified. With them, as in {cmd:gstats winsor}, it lets the new variables
overwrite existing ones (otherwise an existing one exits with r(110)): an
existing target keeps its place in the data and gets the values, storage type,
label and format a new one would, missing outside the sample. A target may not
be another source variable; a source as its own target (for example
{cmd:generate(x) replace} for {cmd:x}) is winsorized in place, but not with
{cmd:if}, {cmd:in} or weights.

{phang}
{opt label} appends the cuts to the target's variable label (its name when
unlabeled), as {cmd:winsor2} writes them: {it:label}{cmd:-Winsor(p}{it:#}{cmd:,p}{it:#}{cmd:)},
or {cmd:-Trim(...)} with {opt trim}, with a 0 before a lower percentile below 1
(for example {cmd:-Winsor(p0.5,p99)}). {cmd:winsor2} adds one more 0 for
each further variable in {it:varlist} ({cmd:p00.5} for the second); {cmd:cwinsor}
writes every label as {cmd:winsor2} writes the first.

{dlgtab:Data recovery}

{pstd}
By default, the command restores the data after an error or interruption.
In-place operations retain a copy of the values being changed until completion.

{phang}
{opt snapshot} uses {help preserve} to copy the whole dataset for recovery.
The default also uses this recovery method when {opt replace} is combined
with {opt generate()}, {opt suffix()}, or {opt prefix()}, so existing output
values, storage types, and metadata survive a failed call.
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

{dlgtab:gtools options}

{pstd}
{cmd:gstats winsor}'s gtools options are accepted so that its scripts run
unchanged. {opt benchmark} and {opt benchmarklevel(#)} display the timing
breakdown of {opt verbose}. {opt compress} and {opt forcestrl} concern how gtools
reads strL variables and have no effect. {opt hashmethod()} ({cmd:0}, {cmd:1},
{cmd:2}, {cmd:default}, {cmd:biject} or {cmd:spooky}; anything else exits with
r(111), as there) and {opt oncollision()} ({cmd:fallback} or {cmd:error}, else
r(198)) choose how gtools hashes the {opt by()} groups, which {cmd:cwinsor} sorts
instead; they are checked and have no effect.


{marker remarks}{...}
{title:Remarks}

{pstd}
Percentiles are computed from nonmissing observations, separately for each
variable and, with {opt by()}, each group. Groups with one or two nonmissing
observations are included. Missing numeric or string group values form groups.

{pstd}
The percentile rule follows {cmd:winsor2}. For sorted nonmissing values x[1],
..., x[N] and p between 0 and 100, let h=N*p/100. If h is an integer, the
cutoff is (x[h]+x[h+1])/2; otherwise it is x[ceil(h)]. The endpoints are the
minimum and maximum. Without {opt by()}, storing the cutoffs in Stata local
macros can round their values.

{pstd}
With weights, the cutoffs are weighted percentiles, as {cmd:gstats winsor}
computes them (and {cmd:_pctile} with weights): with W[i] the cumulative weight
of the sorted values and P = W[N]*p/100, the cutoff is the first x[i] whose
W[i] reaches P, or (x[i]+x[i+1])/2 when W[i] equals P; the endpoints are the
minimum and maximum. Observations with missing or zero weights are left out (as
outside {cmd:if}/{cmd:in}); negative weights exit with r(402), also for
{opt iweight}s, which {cmd:gstats winsor} accepts, and noninteger
{opt fweight}s with r(401). The weighted cutoffs are used as computed, not
rounded through local macros or stored in the source's type with {opt by()}.
{cmd:gstats winsor} sets every value to 0 with weights and {cmd:cuts(0 100)};
{cmd:cwinsor} leaves the values unchanged then. Unlike {cmd:gstats winsor},
the weighted result can replace the variables in place, the default.

{pstd}
New variables from grouped winsorization use Stata's default numeric storage
type. Other outputs retain the source storage type, except when a fractional
cutoff replaces values in an integer variable: its storage type changes from
{cmd:byte} or {cmd:int} to {cmd:float}, or from {cmd:long} to {cmd:double},
following Stata's {cmd:replace}. With {opt by()}, each group's cutoffs are
first stored in the source's type, as {cmd:winsor2} stores them: truncated for
byte, int, and long, and rounded for float. Missing source values do not enter
percentile calculations; generated outputs use system missing for them.


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

{pstd}Weighted percentiles, as {cmd:gstats winsor} takes them, into a named variable:{p_end}
{phang2}{cmd:. cwinsor price [aw=weight], cuts(5 95) generate(price_w5) label}{p_end}

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
{cmd:cwinsor} stores the following in {cmd:r()}, as {cmd:gstats winsor} does:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(J)}}number of groups (1 without {opt by()}){p_end}
{synopt:{cmd:r(minJ)}}observations in the smallest group{p_end}
{synopt:{cmd:r(maxJ)}}observations in the largest group{p_end}
{synopt:{cmd:r(cutlow)}}lower percentile{p_end}
{synopt:{cmd:r(cuthigh)}}upper percentile{p_end}

{pstd}
When the {opt verbose} option is specified, timing information is displayed
but not stored.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
{cmd:cwinsor} implements winsorization similar to Yujun Lian's {cmd:winsor2}.
The original {cmd:winsor2} package is available from SSC:
{browse "https://ideas.repec.org/c/boc/bocode/s457765.html"}.


{title:Also see}

{psee}
Online: {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2} (if installed),
{help summarize}, {help centile}, {help pctile}, {help ctools}
{p_end}

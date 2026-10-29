{smcl}
{* *! version 1.2.0 20261029}{...}
{viewerjumpto "Syntax" "cencode##syntax"}{...}
{viewerjumpto "Description" "cencode##description"}{...}
{viewerjumpto "Options" "cencode##options"}{...}
{viewerjumpto "Remarks" "cencode##remarks"}{...}
{viewerjumpto "Examples" "cencode##examples"}{...}
{viewerjumpto "Stored results" "cencode##results"}{...}
{title:Title}

{phang}
{bf:cencode} {hline 2} C-accelerated parallel string encoding for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cencode}
{varlist}
{ifin}
{cmd:,} {{opth gen:erate(newvarlist)} | {opt replace}} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required (one of)}
{synopt:{opth gen:erate(newvarlist)}}names of new numeric variables{p_end}
{synopt:{opt replace}}replace {it:varlist} with the encoded variables (ctools only){p_end}

{syntab:Options}
{synopt:{opth lab:el(name)}}name for value label (default: same as {it:newvar}){p_end}
{synopt:{opt noext:end}}do not extend existing value label{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cencode} converts string variables to numeric codes with value labels.
It accepts one or more source variables. Specify {opt generate()} with one new
name per source, or {opt replace} to keep the source names.

{pstd}
For a new value label, nonempty strings are sorted in byte order and assigned
codes 1, 2, 3, and so on. Existing labels retain their codes; new strings receive
codes above the largest existing value. If several codes have the same label
text, the smallest is used. Empty strings and observations outside
{cmd:if}/{cmd:in} receive numeric missing values, including with {opt replace}.
Strings such as {cmd:"."} and {cmd:".a"} are encoded as text.


{marker options}{...}
{title:Options}

{dlgtab:Required (one of)}

{phang}
{opth generate(newvarlist)} specifies the names of the new numeric variables to be
created. The number of names must match the number of variables in {it:varlist}.

{phang}
{opt replace} specifies that the variables in {it:varlist} should be replaced with their
encoded numeric versions instead of creating new variables. This is a ctools-specific
extension not available in Stata's built-in {cmd:encode}. When {opt replace}
is specified, each original string variable is dropped and replaced with the
new numeric variable using the same name.

{pstd}
You must specify either {opt generate()} or {opt replace}, but not both.

{dlgtab:Options}

{phang}
{opth label(name)} specifies the name of the value label to be created or
extended. If not specified, the value label is named the same as the new variable.
With multiple source variables, {opt label()} uses and extends the same value
label for each variable in order.

{phang}
{opt noextend} specifies that if the value label already exists, it should not
be extended with additional values. An unlabeled nonempty string causes an
error. By default, new strings are added to the existing label.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
Numeric output uses {cmd:long} storage on the C path. {cmd:encode} handles
{cmd:strL} inputs and existing labels that require long text, extended missing
codes, or codes outside the {cmd:long} range. The value-label limit is 65,536
entries.

{pstd}
Output names must be distinct new variables. If a multi-variable conversion
fails, the original variables and value labels are restored.


{marker examples}{...}
{title:Examples}

{pstd}Encode a string variable and inspect its value label:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cencode make, generate(make_code)}{p_end}
{phang2}{cmd:. label list make_code}{p_end}

{pstd}Use a custom value-label name:{p_end}
{phang2}{cmd:. cencode make, generate(make_num) label(car_makes)}{p_end}

{pstd}Encode selected observations:{p_end}
{phang2}{cmd:. cencode make if foreign == 1, generate(foreign_make)}{p_end}

{pstd}Encode two variables in one call:{p_end}
{phang2}{cmd:. generate str8 origin = cond(foreign, "Foreign", "Domestic")}{p_end}
{phang2}{cmd:. cencode make origin, generate(make_id origin_id) threads(4) verbose}{p_end}

{pstd}Replace the source strings with labeled numeric codes:{p_end}
{phang2}{cmd:. cencode make origin, replace}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cencode} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N_unique)}}number of unique string values encoded (for last variable){p_end}
{synopt:{cmd:r(N_vars)}}number of variables encoded{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] encode}

{psee}
Online: {help encode}, {help decode}, {help cdecode}, {help label}, {help ctools}
{p_end}

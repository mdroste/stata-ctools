{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "cdecode##syntax"}{...}
{viewerjumpto "Description" "cdecode##description"}{...}
{viewerjumpto "Options" "cdecode##options"}{...}
{viewerjumpto "Examples" "cdecode##examples"}{...}
{viewerjumpto "Stored results" "cdecode##results"}{...}
{title:Title}

{phang}
{bf:cdecode} {hline 2} Decode value labels to string variables


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cdecode}
{varlist}
{ifin}
{cmd:,} {{opth gen:erate(newvarlist)} | {opt replace}} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required (one of)}
{synopt:{opth gen:erate(newvarlist)}}create new string variables with decoded labels{p_end}
{synopt:{opt replace}}replace {it:varlist} with the decoded variables (ctools only){p_end}

{syntab:Options}
{synopt:{opth maxl:ength(#)}}maximum string length for new variables{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads for the ctools plugin{p_end}
{synopt:{opt v:erbose}}report how each variable was decoded{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cdecode} extends Stata's {help decode:decode} command to multiple variables
and in-place replacement. It converts numeric variables with value labels
to string variables containing the label text, and produces the same values and
storage types as {cmd:decode}, including literal text and labeled extended
missing values.

{pstd}
Label text is read directly from Stata and written by the ctools plugin in
parallel. Results that the plugin does not produce are decoded by {cmd:decode}
itself: strL results (label text longer than 2,045 bytes), {opt maxlength()}
cuts within text containing multibyte characters, and any call for which the
plugin is unavailable. {cmd:r(N_plugin)} reports how many variables the plugin
decoded.

{pstd}
The source variables must be numeric and have value labels attached. For each
observation, {cmd:cdecode} looks up the numeric value in the value label and
writes the corresponding label text to the new string variable. Unlabeled
values, including unlabeled missing values, produce empty strings.

{pstd}
This is the counterpart to {help cencode:cencode}, which converts string
variables to labeled numeric variables.

{pstd}
Observations outside {cmd:if}/{cmd:in} receive empty strings, including with
{opt replace}. A failed multi-variable conversion restores the original
variables.


{marker options}{...}
{title:Options}

{dlgtab:Required (one of)}

{phang}
{opth generate(newvarlist)} creates new string variables containing the decoded
label text. The number of names must match the number of variables in {it:varlist}.

{phang}
{opt replace} specifies that the variables in {it:varlist} should be replaced with their
decoded string versions instead of creating new variables. The decoded variables keep the original names.

{pstd}
You must specify either {opt generate()} or {opt replace}, but not both.

{dlgtab:Options}

{phang}
{opth maxlength(#)} specifies the maximum string length for the new variables.
By default, the output type follows {cmd:decode}, including strL for
long labels. {cmd:maxlength()} follows {cmd:decode}'s length and truncation rules. Use this option to truncate
labels to a specific length.

{phang}
{opt threads(#)} limits the number of threads the ctools plugin uses. The
argument is an integer from 1 to 256. The default follows the
{help ctools##threads:ctools thread setting}. It has no effect on variables decoded by
{cmd:decode} itself.

{phang}
{opt verbose} reports the number of variables decoded by the ctools plugin and
by {cmd:decode}.


{marker examples}{...}
{title:Examples}

{pstd}Setup: create sample data with value labels{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. gen region = ceil(runiform() * 4)}{p_end}
{phang2}{cmd:. label define region_lbl 1 "North" 2 "South" 3 "East" 4 "West"}{p_end}
{phang2}{cmd:. label values region region_lbl}{p_end}

{pstd}Basic decoding:{p_end}
{phang2}{cmd:. cdecode region, generate(region_str)}{p_end}

{pstd}Limit string length:{p_end}
{phang2}{cmd:. cdecode region, generate(region_short) maxlength(3)}{p_end}

{pstd}Decode with verbose output:{p_end}
{phang2}{cmd:. cdecode region, generate(region_str2) verbose}{p_end}

{pstd}Decode only for certain observations:{p_end}
{phang2}{cmd:. cdecode region if region <= 2, generate(region_ns)}{p_end}

{pstd}Decode two variables in one call:{p_end}
{phang2}{cmd:. clonevar region2 = region}{p_end}
{phang2}{cmd:. cdecode region region2, generate(label1 label2)}{p_end}

{pstd}Replace the source codes with strings:{p_end}
{phang2}{cmd:. cdecode region region2, replace}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cdecode} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N_vars)}}number of variables decoded{p_end}
{synopt:{cmd:r(N_plugin)}}number of those decoded by the ctools plugin{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] decode}

{psee}
Online: {help decode}, {help encode}, {help cencode}, {help label}, {help ctools}
{p_end}

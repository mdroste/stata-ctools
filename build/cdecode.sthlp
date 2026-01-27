{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "cdecode##syntax"}{...}
{viewerjumpto "Description" "cdecode##description"}{...}
{viewerjumpto "Options" "cdecode##options"}{...}
{viewerjumpto "Examples" "cdecode##examples"}{...}
{viewerjumpto "Stored results" "cdecode##results"}{...}
{title:Title}

{phang}
{bf:cdecode} {hline 2} C-accelerated numeric to string decoding for Stata

{pstd}
{cmd:cdecode} is a high-performance drop-in replacement for Stata's {help decode:decode}
command.


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cdecode}
{varname}
{ifin}
{cmd:,} {opth gen:erate(newvar)} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth gen:erate(newvar)}}create new string variable with decoded labels{p_end}

{syntab:Options}
{synopt:{opth maxl:ength(#)}}maximum string length for new variable{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cdecode} is a high-performance drop-in replacement for Stata's {help decode:decode}
command. It converts a numeric variable with value labels to a string variable
containing the label text.

{pstd}
The source variable must be numeric and have a value label attached. For each
observation, {cmd:cdecode} looks up the numeric value in the value label and
writes the corresponding label text to the new string variable.

{pstd}
This is the counterpart to {help cencode:cencode}, which converts string
variables to labeled numeric variables.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth generate(newvar)} creates a new string variable containing the decoded
label text. This option is required.

{dlgtab:Options}

{phang}
{opth maxlength(#)} specifies the maximum string length for the new variable.
By default, {cmd:cdecode} automatically determines the length based on the
longest label in the value label definition. Use this option to truncate
labels to a specific length.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
processing. By default, {cmd:cdecode} uses all available CPU cores.

{phang}
{opt verbose} displays a timing breakdown showing time spent in each phase of
the decoding process.


{marker examples}{...}
{title:Examples}

{pstd}Setup: create sample data with value labels{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000000}{p_end}
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


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cdecode} stores the following in {cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:_cdecode_n_decoded}}number of observations successfully decoded{p_end}
{synopt:{cmd:_cdecode_n_missing}}number of missing values{p_end}
{synopt:{cmd:_cdecode_n_unlabeled}}number of values without labels{p_end}

{pstd}
When the {opt verbose} option is specified, {cmd:cdecode} additionally stores:

{synopt:{cmd:_cdecode_time_parse}}time to parse arguments (seconds){p_end}
{synopt:{cmd:_cdecode_time_decode}}time to decode values (seconds){p_end}
{synopt:{cmd:_cdecode_time_total}}total C plugin time (seconds){p_end}


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

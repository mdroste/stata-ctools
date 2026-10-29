{smcl}
{* *! version 1.2.1 20261029}{...}
{viewerjumpto "Syntax" "cdestring##syntax"}{...}
{viewerjumpto "Description" "cdestring##description"}{...}
{viewerjumpto "Options" "cdestring##options"}{...}
{viewerjumpto "Remarks" "cdestring##remarks"}{...}
{viewerjumpto "Examples" "cdestring##examples"}{...}
{viewerjumpto "Stored results" "cdestring##results"}{...}
{title:Title}

{phang}
{bf:cdestring} {hline 2} C-accelerated parallel string to numeric conversion for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cdestring}
[{varlist}]
{ifin}
{cmd:,} {{opth gen:erate(newvarlist)} | {opt replace}} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required (one of)}
{synopt:{opth gen:erate(newvarlist)}}create new numeric variables{p_end}
{synopt:{opt replace}}replace string variables with numeric versions{p_end}

{syntab:Options}
{synopt:{opth ig:nore(chars)}}remove specified characters before converting{p_end}
{synopt:{opt force}}convert nonnumeric strings to missing{p_end}
{synopt:{opt float}}generate numeric variables as type {cmd:float}{p_end}
{synopt:{opt percent}}convert percentages to fractional form{p_end}
{synopt:{opt dpcomma}}interpret comma as decimal point{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cdestring} converts strings containing numbers to numeric variables.


{pstd}
If {it:varlist} is not specified, {cmd:cdestring} attempts to convert all string
variables in the dataset to numeric.

{pstd}
{cmd:cdestring} creates new numeric variables (when using {opt generate()}) or
replaces the string variables with numeric ones (when using {opt replace}).
Values are parsed in double precision, or float precision with {opt float}.
The output is then compressed to the smallest storage type that retains those
values.

{pstd}
As with {cmd:destring}, empty strings, strings containing only whitespace, and
{cmd:"."} are converted to system missing ({cmd:.}), and {cmd:".a"}, {cmd:".b"},
..., {cmd:".z"} to the corresponding extended missing values. Any other string
that is not a number, including {cmd:"NA"} and {cmd:"NaN"}, is nonnumeric.

{pstd}
If a variable contains any nonnumeric value (after the characters in
{opt ignore()} are removed), {cmd:cdestring} does not generate or replace that
variable unless {opt force} is specified; see {opt force}.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth generate(newvarlist)} creates new numeric variables. The number of names
in {it:newvarlist} must equal the number of variables in {it:varlist}. This
option may not be combined with {opt replace}.

{phang}
{opt replace} replaces string variables with numeric variables of the same name.
This option may not be combined with {opt generate()}.

{dlgtab:Options}

{phang}
{opth ignore(chars)} causes the specified characters to be removed from the
string before attempting numeric conversion. For example, {cmd:ignore("$,")}
removes dollar signs and commas, so "$1,234.56" becomes "1234.56".

{pmore}
Common characters to ignore include currency symbols ({cmd:$}, {cmd:€}, {cmd:£}),
thousands separators ({cmd:,}), and whitespace.

{phang}
{opt force} forces conversion of strings that contain nonnumeric characters not
specified in {opt ignore()}. Such strings are converted to missing values.
Without {opt force}, a variable containing any such string is left unchanged:
no variable is generated for it and it is not replaced. A message such as
{cmd:x: contains nonnumeric characters; no generate} is displayed, the other
variables are still converted, and the return code is 0, as with
{cmd:destring}. Numbers outside the range of a {cmd:double} (for example,
{cmd:"1e400"}), or of a {cmd:float} with {opt float}, count as nonnumeric.

{phang}
{opt float} generates the new numeric variables as type {cmd:float} rather than
double precision before compression. This uses less memory but has less precision.

{phang}
{opt percent} handles percentage values. It removes any percent signs ({cmd:%})
found in the strings, and if any value of a variable contains {cmd:%}, all
values of that variable are divided by 100, as with {cmd:destring}. For
example, "50%" becomes 0.50, and "25" in the same variable becomes 0.25.
{opt percent} implies {cmd:ignore("%")}.

{phang}
{opt dpcomma} specifies that the decimal point is represented by a comma. As
with {cmd:destring}, the first comma is taken as the decimal point, and a value
containing a period is nonnumeric; for example, "3,5" is 3.5, but "3.5" is
nonnumeric. To convert values that use a period as the thousands separator,
such as "1.234,56", also specify {cmd:ignore(".")}.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt verbose} displays progress and elapsed times.


{marker remarks}{...}
{title:Remarks}

{pstd}
Both fixed-width strings ({cmd:str#}) and long strings ({cmd:strL}) are supported.
Only observations selected by {cmd:if}/{cmd:in} are converted. Generated
variables are missing outside that sample. With {opt replace}, unselected
observations also become numeric missing because the source variable changes
from string to numeric. Already numeric variables are skipped.


{pstd}
{bf:Numeric formats recognized:}

{p 8 12 2}- Integers: "123", "-456", "+789"{p_end}
{p 8 12 2}- Decimals: "123.456", "-0.123", ".5"{p_end}
{p 8 12 2}- Scientific notation: "1.23e4", "1.23E-4", "1e10" (the exponent needs a digit: "1e+" is nonnumeric){p_end}
{p 8 12 2}- Decimal comma (with {opt dpcomma}): "1234,56"; with {cmd:ignore(".")} also "1.234,56"{p_end}

{pstd}
{bf:Missing value representations:}

{p 8 12 2}- Empty string or whitespace only: {cmd:.}{p_end}
{p 8 12 2}- Single period: "." becomes {cmd:.}{p_end}
{p 8 12 2}- ".a" to ".z" become {cmd:.a} to {cmd:.z} (with {opt percent}, they become {cmd:.} when the variable is divided by 100){p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup: create sample data with string numbers{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. gen str10 price_str = string(runiform()*1000, "%9.2f")}{p_end}
{phang2}{cmd:. gen str10 qty_str = string(round(runiform()*100))}{p_end}

{pstd}Basic conversion with generate:{p_end}
{phang2}{cmd:. cdestring price_str qty_str, generate(price qty)}{p_end}

{pstd}Handle currency formatting:{p_end}
{phang2}{cmd:. gen str15 cost = "$" + string(runiform()*10000, "%10.2fc")}{p_end}
{phang2}{cmd:. cdestring cost, replace ignore("$,")}{p_end}

{pstd}Handle European number format (period as thousands separator):{p_end}
{phang2}{cmd:. gen str12 euro_price = "1.234,56"}{p_end}
{phang2}{cmd:. cdestring euro_price, generate(euro_num) dpcomma ignore(".")}{p_end}

{pstd}Convert percentages:{p_end}
{phang2}{cmd:. gen str6 pct = string(runiform()*100, "%5.1f") + "%"}{p_end}
{phang2}{cmd:. cdestring pct, generate(pct_frac) percent}{p_end}

{pstd}Force conversion of mixed data:{p_end}
{phang2}{cmd:. replace price_str = "N/A" in 1/10}{p_end}
{phang2}{cmd:. cdestring price_str, generate(price_forced) force}{p_end}

{pstd}Convert with verbose timing output:{p_end}
{phang2}{cmd:. cdestring qty_str, generate(qty_verbose) verbose}{p_end}

{pstd}Limit to 4 threads:{p_end}
{phang2}{cmd:. cdestring qty_str, generate(qty_parallel) threads(4)}{p_end}

{pstd}Convert only for certain observations:{p_end}
{phang2}{cmd:. cdestring price_str if qty > 50, generate(price_filtered) force}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cdestring} stores the following in global scalars (accessible via
{cmd:scalar(_cdestring_...)}):

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:_cdestring_n_converted}}number of values that were numbers or missing-value strings{p_end}
{synopt:{cmd:_cdestring_n_failed}}number of nonnumeric values, summed over variables (with or without {opt force}){p_end}

{pstd}
When the {opt verbose} option is specified, {cmd:cdestring} additionally stores:

{synopt:{cmd:_cdestring_time_parse}}time to parse arguments (seconds){p_end}
{synopt:{cmd:_cdestring_time_convert}}time to convert strings (seconds){p_end}
{synopt:{cmd:_cdestring_time_total}}total C plugin time (seconds){p_end}
{synopt:{cmd:_cdestring_openmp_enabled}}1 if OpenMP is enabled, 0 otherwise{p_end}
{synopt:{cmd:_cdestring_threads_max}}maximum available threads{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] destring}

{psee}
Online: {help destring}, {help tostring}, {help encode}, {help real()}, {help ctools}
{p_end}

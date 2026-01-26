{smcl}
{* *! version 0.9.0 26Jan2026}{...}
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
{cmd:cdestring} is a high-performance drop-in replacement for Stata's {help destring:destring}
command. It converts string variables containing numeric values to numeric variables,
using parallel algorithms for significantly improved performance on large datasets.

{pstd}
If {it:varlist} is not specified, {cmd:cdestring} attempts to convert all string
variables in the dataset to numeric.

{pstd}
{cmd:cdestring} creates new numeric variables (when using {opt generate()}) or
replaces the string variables with numeric ones (when using {opt replace}).
By default, variables are created as type {cmd:double}; use {opt float} for
single-precision floating point.

{pstd}
Empty strings and strings containing only whitespace are converted to missing
values. The standard Stata missing value representations ({cmd:.}, {cmd:NA},
{cmd:NaN}) are also converted to missing values.


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
Without {opt force}, the conversion proceeds but a warning is issued.

{phang}
{opt float} generates the new numeric variables as type {cmd:float} rather than
the default {cmd:double}. This uses less memory but has less precision.

{phang}
{opt percent} handles percentage values. It removes any percent signs ({cmd:%})
found in the strings and divides the values by 100 to convert them to fractional
form. For example, "50%" becomes 0.50.

{phang}
{opt dpcomma} specifies that the decimal point is represented by a comma and
the grouping separator is a period (European format). For example, "1.234,56"
is interpreted as one thousand two hundred thirty-four and 56/100.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
conversion. By default, {cmd:cdestring} uses all available CPU cores as reported
by OpenMP. Use this option to limit parallelism.

{phang}
{opt verbose} displays a timing breakdown showing time spent in each phase of
the conversion process.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cdestring} uses a high-performance C plugin that parallelizes the conversion
process. The algorithm works as follows:

{p 8 12 2}1. {bf:Parallel processing:} Each observation is processed independently
across multiple CPU cores using OpenMP.{p_end}

{p 8 12 2}2. {bf:Character stripping:} If {opt ignore()} is specified, the
specified characters are removed using an O(1) lookup table.{p_end}

{p 8 12 2}3. {bf:Fast parsing:} A highly optimized numeric parser converts
strings to doubles, handling integers, decimals, and scientific notation.{p_end}

{p 8 12 2}4. {bf:Parallel storage:} Results are written back to Stata in
parallel when the dataset is large enough.{p_end}

{pstd}
For large datasets (millions of observations), {cmd:cdestring} can be significantly
faster than the native {cmd:destring} command.

{pstd}
{cmd:cdestring} supports both regular string variables ({cmd:str#}) and long strings
({cmd:strL}).

{pstd}
{bf:Numeric formats recognized:}

{p 8 12 2}- Integers: "123", "-456", "+789"{p_end}
{p 8 12 2}- Decimals: "123.456", "-0.123", ".5"{p_end}
{p 8 12 2}- Scientific notation: "1.23e4", "1.23E-4", "1e10"{p_end}
{p 8 12 2}- European format (with {opt dpcomma}): "1.234,56"{p_end}

{pstd}
{bf:Missing value representations:}

{p 8 12 2}- Empty string or whitespace only{p_end}
{p 8 12 2}- Single period: "."{p_end}
{p 8 12 2}- "NA", "na"{p_end}
{p 8 12 2}- "NaN", "nan"{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup: create sample data with string numbers{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000000}{p_end}
{phang2}{cmd:. gen str10 price_str = string(runiform()*1000, "%9.2f")}{p_end}
{phang2}{cmd:. gen str10 qty_str = string(round(runiform()*100))}{p_end}

{pstd}Basic conversion with generate:{p_end}
{phang2}{cmd:. cdestring price_str qty_str, generate(price qty)}{p_end}

{pstd}Replace the original string variables:{p_end}
{phang2}{cmd:. cdestring price_str qty_str, replace}{p_end}

{pstd}Handle currency formatting:{p_end}
{phang2}{cmd:. gen str15 cost = "$" + string(runiform()*10000, "%10.2fc")}{p_end}
{phang2}{cmd:. cdestring cost, replace ignore("$,")}{p_end}

{pstd}Handle European number format:{p_end}
{phang2}{cmd:. gen str12 euro_price = "1.234,56"}{p_end}
{phang2}{cmd:. cdestring euro_price, generate(euro_num) dpcomma}{p_end}

{pstd}Convert percentages:{p_end}
{phang2}{cmd:. gen str6 pct = string(runiform()*100, "%5.1f") + "%"}{p_end}
{phang2}{cmd:. cdestring pct, generate(pct_frac) percent}{p_end}

{pstd}Force conversion of mixed data:{p_end}
{phang2}{cmd:. replace price_str = "N/A" in 1/10}{p_end}
{phang2}{cmd:. cdestring price_str, replace force}{p_end}

{pstd}Convert with verbose timing output:{p_end}
{phang2}{cmd:. cdestring price_str, replace verbose}{p_end}

{pstd}Limit to 4 threads:{p_end}
{phang2}{cmd:. cdestring price_str qty_str, replace threads(4)}{p_end}

{pstd}Convert only for certain observations:{p_end}
{phang2}{cmd:. cdestring price_str if qty > 50, generate(price_filtered)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cdestring} stores the following in {cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:_cdestring_n_converted}}number of observations successfully converted{p_end}
{synopt:{cmd:_cdestring_n_failed}}number of observations with nonnumeric values{p_end}

{pstd}
When the {opt verbose} option is specified, {cmd:cdestring} additionally stores:

{synopt:{cmd:_cdestring_time_parse}}time to parse arguments (seconds){p_end}
{synopt:{cmd:_cdestring_time_convert}}time to convert strings (seconds){p_end}
{synopt:{cmd:_cdestring_time_total}}total C plugin time (seconds){p_end}
{synopt:{cmd:_cdestring_openmp_enabled}}1 if OpenMP is enabled, 0 otherwise{p_end}
{synopt:{cmd:_cdestring_threads_max}}maximum available threads{p_end}


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help destring}, {help tostring}, {help encode}, {help real()}, {help ctools}
{p_end}

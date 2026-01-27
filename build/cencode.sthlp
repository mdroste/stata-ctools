{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "cencode##syntax"}{...}
{viewerjumpto "Description" "cencode##description"}{...}
{viewerjumpto "Options" "cencode##options"}{...}
{viewerjumpto "Remarks" "cencode##remarks"}{...}
{viewerjumpto "Examples" "cencode##examples"}{...}
{viewerjumpto "Stored results" "cencode##results"}{...}
{title:Title}

{phang}
{bf:cencode} {hline 2} C-accelerated parallel string encoding for Stata

{pstd}
{cmd:cencode} is a high-performance drop-in replacement for Stata's {help encode:encode}
command.


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cencode}
{varname}
{ifin}
{cmd:,} {opth gen:erate(newvar)} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth gen:erate(newvar)}}name of new numeric variable{p_end}

{syntab:Options}
{synopt:{opth lab:el(name)}}name for value label (default: same as {it:newvar}){p_end}
{synopt:{opt noext:end}}do not extend existing value label{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cencode} is a high-performance drop-in replacement for Stata's {help encode:encode}
command. It converts a string variable to a numeric variable with a corresponding
value label, using parallel algorithms for significantly improved performance
on large datasets.

{pstd}
{cmd:cencode} creates a new variable named {it:newvar} of type {cmd:long} that contains
an integer code for each unique string value in {varname}. Integer codes are assigned
in alphabetical order of the string values (1, 2, 3, ...). A value label with the
same name as {it:newvar} (or specified by {opt label()}) is created that maps
integer codes back to the original string values.

{pstd}
Missing values (empty strings and extended missing values) in {varname} become
system missing values in {it:newvar}.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth generate(newvar)} specifies the name of the new numeric variable to be
created. This option is required.

{dlgtab:Options}

{phang}
{opth label(name)} specifies the name of the value label to be created or
extended. If not specified, the value label is named the same as {it:newvar}.

{phang}
{opt noextend} specifies that if the value label already exists, it should not
be extended with additional values. By default, {cmd:cencode} adds new string
values to existing labels.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:cencode} uses all available CPU cores as reported by
OpenMP. Use this option to limit parallelism.

{phang}
{opt verbose} displays a timing breakdown showing time spent in each phase of
the encoding: loading strings, collecting unique values, sorting, encoding,
and storing value labels.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cencode} uses a high-performance C plugin that parallelizes the encoding
process. The algorithm works as follows:

{p 8 12 2}1. {bf:Parallel string loading:} String values are loaded from Stata
into C memory using parallel I/O when the dataset is large enough.{p_end}

{p 8 12 2}2. {bf:Unique value collection:} A hash table is used to efficiently
identify all unique string values in O(n) time.{p_end}

{p 8 12 2}3. {bf:Alphabetical sorting:} Unique values are sorted alphabetically
to ensure consistent integer code assignment across runs.{p_end}

{p 8 12 2}4. {bf:Parallel encoding:} Integer codes are written to the new
variable using parallel writes when the dataset is large enough.{p_end}

{p 8 12 2}5. {bf:Value label creation:} The string-to-integer mapping is passed
to Stata to create the value label definition.{p_end}

{pstd}
For large datasets (millions of observations), {cmd:cencode} can be significantly
faster than the native {cmd:encode} command.

{pstd}
{cmd:cencode} supports both regular string variables ({cmd:str#}) and long strings
({cmd:strL}). The maximum number of unique values is 65,536.


{marker examples}{...}
{title:Examples}

{pstd}Basic encoding of a string variable:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cencode make, generate(make_code)}{p_end}

{pstd}View the created variable and its labels:{p_end}
{phang2}{cmd:. tabulate make_code}{p_end}
{phang2}{cmd:. label list make_code}{p_end}

{pstd}Encode with a custom label name:{p_end}
{phang2}{cmd:. cencode make, generate(make_num) label(car_makes)}{p_end}

{pstd}Encode with verbose timing output:{p_end}
{phang2}{cmd:. cencode make, generate(make_code) verbose}{p_end}

{pstd}Limit to 4 threads:{p_end}
{phang2}{cmd:. cencode make, generate(make_code) threads(4)}{p_end}

{pstd}Encode only for certain observations:{p_end}
{phang2}{cmd:. cencode make if foreign == 1, generate(foreign_make)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cencode} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:_cencode_n_unique}}number of unique string values encoded{p_end}

{pstd}
When the {opt verbose} option is specified, {cmd:cencode} additionally stores:

{synopt:{cmd:_cencode_time_load}}time to load strings from Stata (seconds){p_end}
{synopt:{cmd:_cencode_time_collect}}time to collect unique values (seconds){p_end}
{synopt:{cmd:_cencode_time_sort}}time to sort unique values (seconds){p_end}
{synopt:{cmd:_cencode_time_encode}}time to encode values (seconds){p_end}
{synopt:{cmd:_cencode_time_labels}}time to store label mapping (seconds){p_end}
{synopt:{cmd:_cencode_time_total}}total C plugin time (seconds){p_end}
{synopt:{cmd:_cencode_openmp_enabled}}1 if OpenMP is enabled, 0 otherwise{p_end}
{synopt:{cmd:_cencode_threads_max}}maximum available threads{p_end}


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

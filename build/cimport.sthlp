{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "cimport##syntax"}{...}
{viewerjumpto "Description" "cimport##description"}{...}
{viewerjumpto "Options" "cimport##options"}{...}
{viewerjumpto "Remarks" "cimport##remarks"}{...}
{viewerjumpto "Examples" "cimport##examples"}{...}
{viewerjumpto "Stored results" "cimport##results"}{...}
{title:Title}

{phang}
{bf:cimport delimited} {hline 2} C-accelerated CSV/delimited text import


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cimport}
{cmd:delimited}
{cmd:using}
{it:filename}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt d:elimiters(chars)}}specify field delimiter; default is comma{p_end}

{syntab:Variable names}
{synopt:{opt varn:ames(rule)}}rule for reading variable names; {opt 1} or {opt nonames}{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}

{syntab:Parsing}
{synopt:{opt bindq:uotes(option)}}quote binding rule; {opt strict} or {opt loose}{p_end}
{synopt:{opt stripq:uotes}}remove surrounding quotes from string values{p_end}
{synopt:{opt enc:oding(encoding)}}file encoding; currently only UTF-8 supported{p_end}
{synopt:{opt rowr:ange([start][:end])}}range of rows to import{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cimport delimited} is a high-performance replacement for
{help import delimited:import delimited} that uses a C plugin with
multi-threaded parallel parsing.

{pstd}
The command reads delimited text files (CSV, TSV, etc.) and loads them into
Stata, automatically inferring variable types and handling quoted fields.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt clear} specifies that it is okay to replace the data in memory, even
though the current data have not been saved to disk.

{phang}
{opt delimiters(chars)} specifies the delimiter used in the file. The default
is comma ({cmd:,}). Use {cmd:delimiters(tab)} or {cmd:delimiters(\t)} for
tab-delimited files.

{dlgtab:Variable names}

{phang}
{opt varnames(rule)} specifies how variable names are determined.
{opt varnames(1)} treats the first row as variable names (the default).
{opt varnames(nonames)} treats the first row as data and generates default
variable names (v1, v2, ...).

{phang}
{opt case(option)} specifies the case of variable names. {opt preserve}
(the default) keeps the original case. {opt lower} converts to lowercase.
{opt upper} converts to uppercase.

{dlgtab:Parsing}

{phang}
{opt bindquotes(option)} specifies how quoted fields are handled.
{opt strict} (the default) requires fields to be properly quoted.
{opt loose} is more permissive with quote handling.

{phang}
{opt stripquotes} removes surrounding quotation marks from string values
after parsing.

{phang}
{opt encoding(encoding)} specifies the file encoding. Currently only UTF-8
is supported.

{phang}
{opt rowrange([start][:end])} specifies a range of rows to import.
Use {opt rowrange(100:200)} to import rows 100-200, or {opt rowrange(100:)}
to import from row 100 to the end.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information including timing
breakdown and throughput in MB/s.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cimport} uses a three-phase approach:

{p 8 12 2}1. {bf:Scan:} Parse the file to determine column types and widths{p_end}
{p 8 12 2}2. {bf:Create:} Create variables with appropriate Stata types{p_end}
{p 8 12 2}3. {bf:Load:} Load data into variables using parallel processing{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Import a CSV file:{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear}{p_end}

{pstd}Import a tab-delimited file:{p_end}
{phang2}{cmd:. cimport delimited using data.tsv, clear delimiters(tab)}{p_end}

{pstd}Import with verbose output and lowercase variable names:{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear case(lower) verbose}{p_end}

{pstd}Import only rows 1000-2000:{p_end}
{phang2}{cmd:. cimport delimited using bigdata.csv, clear rowrange(1000:2000)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cimport delimited} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations imported{p_end}
{synopt:{cmd:r(k)}}number of variables created{p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(filename)}}name of the imported file{p_end}


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help import delimited}, {help ctools}, {help cexport}
{p_end}

{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "cimport##syntax"}{...}
{viewerjumpto "Description" "cimport##description"}{...}
{viewerjumpto "Options" "cimport##options"}{...}
{viewerjumpto "Remarks" "cimport##remarks"}{...}
{viewerjumpto "Examples" "cimport##examples"}{...}
{viewerjumpto "Stored results" "cimport##results"}{...}
{title:Title}

{phang}
{bf:cimport delimited} {hline 2} C-accelerated CSV/delimited text import

{pstd}
{cmd:cimport delimited} is a high-performance replacement for Stata's
{help import delimited:import delimited} command.


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cimport}
{cmd:delimited}
{cmd:using}
{it:filename}
[{cmd:,} {it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt d:elimiters(chars)}}specify field delimiter; default is comma{p_end}

{syntab:Variable names}
{synopt:{opt varn:ames(rule)}}rule for reading variable names; {opt 1} or {opt nonames}{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}

{syntab:Variable types}
{synopt:{opt asfloat}}import all numeric variables as float{p_end}
{synopt:{opt asdouble}}import all numeric variables as double{p_end}
{synopt:{opt numeric:cols(numlist)}}force specified columns to be numeric{p_end}
{synopt:{opt string:cols(numlist)}}force specified columns to be string{p_end}

{syntab:Parsing}
{synopt:{opt bindq:uotes(option)}}quote binding rule; {opt strict} or {opt loose}{p_end}
{synopt:{opt stripq:uotes}}remove surrounding quotes from string values{p_end}
{synopt:{opt enc:oding(encoding)}}file encoding; currently only UTF-8 supported{p_end}
{synopt:{opt rowr:ange([start][:end])}}range of rows to import{p_end}
{synopt:{opt colr:ange([start][:end])}}range of columns to import{p_end}
{synopt:{opt empty:lines(option)}}empty line handling; {opt skip} or {opt fill}{p_end}

{syntab:Number formats}
{synopt:{opt decimals:eparator(char)}}decimal point character; default is period{p_end}
{synopt:{opt groups:eparator(char)}}thousands grouping character; default is none{p_end}
{synopt:{opt loc:ale(name)}}locale for number parsing (e.g., de_DE, fr_FR){p_end}
{synopt:{opt parsel:ocale}}enable locale-aware number parsing{p_end}
{synopt:{opt maxquoted:rows(#)}}max rows to scan for quote inference; default 20{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synopt:{opt thr:eads(#)}}number of threads to use; default is auto-detect{p_end}
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

{dlgtab:Variable types}

{phang}
{opt asfloat} imports all numeric variables as Stata {help data types:float}
type, regardless of the detected optimal type. This uses less memory than
double but has less precision. Cannot be combined with {opt asdouble}.

{phang}
{opt asdouble} imports all numeric variables as Stata {help data types:double}
type, regardless of whether a smaller type would suffice. This ensures maximum
precision. Cannot be combined with {opt asfloat}.

{phang}
{opt numericcols(numlist)} forces the specified columns to be imported as
numeric. Column numbers are 1-based. Values that cannot be parsed as numbers
become missing. This overrides automatic type detection for these columns.

{phang}
{opt stringcols(numlist)} forces the specified columns to be imported as
string, even if they contain only numeric values. Column numbers are 1-based.
This overrides automatic type detection for these columns.

{dlgtab:Parsing}

{phang}
{opt bindquotes(option)} specifies how quoted fields are handled.
{opt loose} (the default) treats each line as a row, ignoring quotes.
{opt strict} respects quotes so that quoted fields can span multiple lines.

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

{phang}
{opt colrange([start][:end])} specifies a range of columns to import.
Use {opt colrange(2:5)} to import columns 2-5, or {opt colrange(3:)}
to import from column 3 to the last column.

{phang}
{opt emptylines(option)} specifies how empty lines in the file are handled.
{opt skip} (the default) ignores empty lines. {opt fill} includes empty
lines as observations with all missing values.

{dlgtab:Number formats}

{phang}
{opt decimalseparator(char)} specifies the character used as the decimal
point in numeric values. The default is period ({cmd:.}). For European-format
files that use comma as the decimal separator, specify {opt decimalseparator(,)}.

{phang}
{opt groupseparator(char)} specifies the character used as a thousands
grouping separator in numeric values. The default is no grouping separator.
For files with numbers like "1,234,567" use {opt groupseparator(,)}, or for
European formats like "1.234.567" use {opt groupseparator(.)}.

{phang}
{opt locale(name)} specifies the locale to use for parsing numeric values.
Common locales include {opt de_DE} (German), {opt fr_FR} (French),
{opt en_US} (US English). The locale determines the default decimal and
grouping separators when {opt parselocale} is also specified.

{phang}
{opt parselocale} enables locale-aware parsing of numeric values. When
specified, numbers are parsed according to the locale's conventions. If
{opt locale()} is also specified, that locale is used; otherwise the
system default is assumed. This option automatically sets the decimal and
grouping separators based on the locale (e.g., German locale uses comma
as decimal separator and period as grouping separator).

{phang}
{opt maxquotedrows(#)} specifies the maximum number of rows to scan when
inferring whether fields contain quoted strings. The default is 20. Increase
this value if your file has quoted fields that don't appear until later rows.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information including timing
breakdown and throughput in MB/s.

{phang}
{opt threads(#)} specifies the number of threads to use for parallel
parsing. The default (0) auto-detects based on available CPU cores.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cimport} uses a three-phase approach:

{p 8 12 2}1. {bf:Scan:} Parse the file to determine column types and widths{p_end}
{p 8 12 2}2. {bf:Create:} Create variables with appropriate Stata types{p_end}
{p 8 12 2}3. {bf:Load:} Load data into variables using parallel processing{p_end}

{pstd}
{bf:European number formats:} When importing files that use European number
conventions (comma as decimal separator, period as thousands separator),
use both {opt decimalseparator(,)} and {opt groupseparator(.)} together.


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

{pstd}Import European-format file (comma decimal, period grouping):{p_end}
{phang2}{cmd:. cimport delimited using european.csv, clear decimalseparator(,) groupseparator(.)}{p_end}

{pstd}Import German-format file using locale:{p_end}
{phang2}{cmd:. cimport delimited using german.csv, clear locale(de_DE) parselocale}{p_end}

{pstd}Import only columns 2-4:{p_end}
{phang2}{cmd:. cimport delimited using wide.csv, clear colrange(2:4)}{p_end}

{pstd}Force all numerics to double precision:{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear asdouble}{p_end}

{pstd}Force column 3 to be string (e.g., ZIP codes with leading zeros):{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear stringcols(3)}{p_end}

{pstd}Force columns 2 and 4 to be numeric:{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear numericcols(2 4)}{p_end}

{pstd}Include empty lines as missing observations:{p_end}
{phang2}{cmd:. cimport delimited using data.csv, clear emptylines(fill)}{p_end}


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
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] import delimited}

{psee}
Online: {help import delimited}, {help insheet}, {help cexport}, {help ctools}
{p_end}

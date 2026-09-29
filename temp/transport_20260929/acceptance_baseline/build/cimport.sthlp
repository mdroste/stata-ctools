{smcl}
{* *! version 1.0.1 07Feb2026}{...}
{viewerjumpto "Syntax" "cimport##syntax"}{...}
{viewerjumpto "Description" "cimport##description"}{...}
{viewerjumpto "Options" "cimport##options"}{...}
{viewerjumpto "Remarks" "cimport##remarks"}{...}
{viewerjumpto "Examples" "cimport##examples"}{...}
{viewerjumpto "Stored results" "cimport##results"}{...}
{title:Title}

{phang}
{bf:cimport} {hline 2} C-accelerated text, Excel, and statistical file import


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cimport}
{cmd:delimited}
[{cmd:using}]
{it:filename}
[{cmd:,} {it:options}]

{p 8 17 2}
{cmdab:cimport}
{cmd:delimited}
{it:extvarlist}
{cmd:using}
{it:filename}
[{cmd:,} {it:options}]

{pstd}
As with {cmd:import delimited}, {it:extvarlist} renames the first imported
variables, in order and exactly as typed, after the import; variable names are
still read from the file as usual. Specifying more names than imported
variables is an error.

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt d:elimiters(chars)}}specify field delimiter; default is automatic detection{p_end}

{syntab:Variable names}
{synopt:{opt varn:ames(rule)}}rule for reading variable names; {opt 1} or {opt nonames}{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}

{syntab:Variable types}
{synopt:{opt asfloat}}import all numeric variables as float{p_end}
{synopt:{opt asdouble}}import all numeric variables as double{p_end}
{synopt:{opt numeric:cols(numlist)}}force specified columns to be numeric{p_end}
{synopt:{opt string:cols(numlist)}}force specified columns to be string{p_end}

{syntab:Parsing}
{synopt:{opt bindq:uotes(option)}}quote binding rule; {opt strict}, {opt loose}, or {opt nobind}{p_end}
{synopt:{opt stripq:uotes(default|yes|no)}}control removal of surrounding quotes{p_end}
{synopt:{opt enc:oding(encoding)}}file encoding; default is automatic detection; UTF-8/16/32 and named system encodings{p_end}
{synopt:{opt rowr:ange([start][:end])}}range of rows to import{p_end}
{synopt:{opt colr:ange([start][:end])}}range of columns to import{p_end}
{synopt:{opt empty:lines(option)}}empty line handling; {opt skip} or {opt fill}{p_end}

{syntab:Number formats}
{synopt:{opt decimals:eparator(char)}}decimal point character; default is period{p_end}
{synopt:{opt groups:eparator(char)}}thousands grouping character; default is none{p_end}
{synopt:{opt loc:ale(name)}}locale for number parsing (e.g., de_DE, fr_FR){p_end}
{synopt:{opt parsel:ocale(name)}}locale for numeric parsing{p_end}
{synopt:{opt maxquoted:rows(#)}}limit on lines in a strict quoted field{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{pstd}Import Excel (.xls or .xlsx) file:

{p 8 17 2}
{cmdab:cimport}
{cmd:excel}
[{cmd:using}]
{it:filename}
[{cmd:,} {it:excel_options}]

{synoptset 32 tabbed}{...}
{synopthdr:excel options}
{synoptline}
{syntab:Main}
{synopt:{opt sheet(name)}}worksheet name to import{p_end}
{synopt:{opt cellr:ange(range)}}cell range to import (e.g., A1:F100){p_end}
{synopt:{opt first:row}}treat first row as variable names{p_end}
{synopt:{opt alls:tring}}import all columns as string{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt desc:ribe}}list worksheet names and ranges without loading data{p_end}
{synopt:{opt v:erbose}}display progress information{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cimport} reads delimited text, Excel (.xls and .xlsx), dBase, SAS7BDAT
with optional SAS7BCAT catalogs, SAS XPORT5/8, SPSS SAV/ZSAV, and ESRI shapefiles.
File parsing and I/O run in C. FRED and Windows Haver providers are excluded from this parity scope.

{phang2}
{cmd:cimport delimited} imports delimited text files (CSV, TSV, etc.).
It is a high-performance replacement for {help import delimited:import delimited}.

{phang2}
{cmd:cimport excel} imports Excel (.xls and .xlsx) files. It is a high-performance
replacement for {help import excel:import excel}.

{pstd}
Other file formats use {cmd:cimport dbase}, {cmd:cimport sas},
{cmd:cimport sasxport5}, {cmd:cimport sasxport8}, {cmd:cimport spss}, and
{cmd:cimport shp}, followed by {cmd:using filename}. SAS/SPSS also accept a
variable list and if/in selection. SAS accepts {cmd:bcat(filename)} and
{cmd:encoding()}; SPSS accepts {cmd:encoding()} and {cmd:zsav}; XPORT5 accepts
{cmd:member()}, {cmd:novallabels}, and {cmd:describe}.

{pstd}
Excel column lists may use {cmd:identifier=A note=C}; {cmd:describe} returns
{cmd:r(N_worksheet)}, {cmd:r(worksheet_1)}, {cmd:r(range_1)}, and corresponding
results for further sheets. {cmd:allstring} accepts an optional numeric format;
locale options are accepted. See the audit for Unicode locale validation limits.

{pstd}
Delimited/Excel/statistical readers use strL for text longer than 2,045 bytes;
SHP binary headers also use strL. Native invalid long-string display formats are
replaced with valid {cmd:%9s}; see the parity audit for this metadata difference.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt clear} specifies that it is okay to replace the data in memory, even
though the current data have not been saved to disk.

{phang}
{opt delimiters(chars)} specifies the delimiter used in the file. The default is automatic detection; use {cmd:delimiters(",")} to force comma. Use {cmd:delimiters(tab)} or {cmd:delimiters(\t)} for
tab-delimited files. As in {cmd:import delimited}, automatic detection picks
the most frequent of tab, comma, semicolon, colon, and pipe outside quotes in
the first 50 nonempty lines, including the variable-names line; ties favor
comma, then tab, pipe, colon, and semicolon.

{dlgtab:Variable names}

{phang}
{opt varnames(rule)} specifies how variable names are determined.
By default, headers are inferred. {opt varnames(1)} explicitly treats the first row as variable names.
{opt varnames(nonames)} treats the first row as data and generates default
variable names (v1, v2, ...). You can also specify {opt varnames(}{it:N}{opt )}
to use row {it:N} as variable names (e.g., {opt varnames(3)} uses the third row).

{phang}
{opt case(option)} specifies the case of variable names. Delimited import defaults
to {opt lower}; Excel defaults to {opt preserve}. {opt lower} converts to lowercase.
{opt upper} converts to uppercase.

{dlgtab:Variable types}

{phang}
{opt asfloat} uses Stata {help data types:float} for noninteger numeric
columns. Integer columns retain their inferred byte, int, or long storage. This uses less memory than
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
{opt stripquotes(default)} removes enclosing quotes when they bind a field.
{opt stripquotes(yes)} removes all quotation marks; {opt stripquotes(no)} retains them.

{phang}
{opt encoding(encoding)} overrides automatic encoding detection. Supported encodings include UTF-8,
UTF-16LE/BE, UTF-32LE/BE, ASCII, Latin-1/9, Windows code pages and Mac Roman.
Named legacy encodings use the system C conversion API.

{phang}
{opt rowrange([start][:end])} specifies a range of rows to import. As in
{cmd:import delimited}, rows are lines of the file: the variable-names line,
blank lines, and lines inside quoted fields all count. Use
{opt rowrange(100:200)} to import lines 100-200, and {opt rowrange(100:)} or
{opt rowrange(100)} to import from line 100 to the end.

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
{opt parselocale(name)} selects a supported numeric locale, for example
{cmd:parselocale(de_DE)}. Explicit decimal/group separator options take precedence.
Delimited import uses {opt parselocale(name)}; {opt locale(name)} is not a native delimited option.

{phang}
{opt maxquotedrows(#)} limits physical rows within a strict quoted field and
accepts an integer or {cmd:unlimited}; exceeding it returns 5101.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:cimport} uses all available CPU cores.


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

{pstd}Setup: create a CSV file to import:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. export delimited using auto.csv, replace}{p_end}

{pstd}Import the CSV file:{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear}{p_end}

{pstd}Import with verbose output and lowercase variable names:{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear case(lower) verbose}{p_end}

{pstd}Import only the first 20 lines (the variable names and 19 observations):{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear rowrange(1:20)}{p_end}

{pstd}Import only columns 1-3:{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear colrange(1:3)}{p_end}

{pstd}Force all numerics to double precision:{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear asdouble}{p_end}

{pstd}Force column 1 to be string (make is already string, but demonstrates syntax):{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear stringcols(1)}{p_end}

{pstd}Create and import a tab-delimited file:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. export delimited using auto.tsv, delimiter(tab) replace}{p_end}
{phang2}{cmd:. cimport delimited using auto.tsv, clear delimiters(tab)}{p_end}

{pstd}Import with multi-threading:{p_end}
{phang2}{cmd:. cimport delimited using auto.csv, clear threads(4)}{p_end}

{pstd}{ul:Excel import examples}

{pstd}Import an Excel file:{p_end}
{phang2}{cmd:. cimport excel using mydata.xlsx, clear}{p_end}

{pstd}Import a specific sheet:{p_end}
{phang2}{cmd:. cimport excel using mydata.xlsx, clear sheet("Sheet2")}{p_end}

{pstd}Import with first row as variable names:{p_end}
{phang2}{cmd:. cimport excel using mydata.xlsx, clear firstrow}{p_end}

{pstd}Import all columns as string:{p_end}
{phang2}{cmd:. cimport excel using mydata.xlsx, clear allstring}{p_end}


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
{synopt:{cmd:r(encoding)}}source encoding selected or supplied{p_end}
{synopt:{cmd:r(delimiters)}}field delimiters used; tabs reported as {cmd:\t}{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] import delimited}, {bf:[D] import excel}

{psee}
Online: {help import delimited}, {help import excel}, {help insheet}, {help cexport}, {help ctools}
{p_end}

{title:Validation and failure behavior}

{pstd}
Excel import honors the 1900 and 1904 date systems. As in native import, a
numeric column with any date or time cell becomes a Stata date when all such
cells are dates only and a datetime (milliseconds) otherwise; every numeric
cell in the column, including General-format serials, converts to that unit.
The column keeps the cells' display format only when all its numeric cells
share one; otherwise it uses %td or %tc. In the 1900 system, serial 60 maps to
28 February 1900, following native import. Arbitrary workbook display formats
are not fully mapped to Stata formats.
{p_end}

{pstd}
TRUE/FALSE cells import as 1/0 (a column of them displays as %1.0f), error
cells as missing, and text formula results as strings. XLSX worksheets are
located through the workbook relationships, and the whole worksheet is read
whatever its {it:dimension} record says. Chart, dialog and macro sheets cannot
be imported; {cmd:describe} lists them without a range, as it does for empty
worksheets.
{p_end}

{pstd}
Full feature parity with native import/export is not established. See
{browse "https://github.com/mdroste/stata-ctools/blob/main/docs/IO_PARITY.md":the I/O parity audit}
for unsupported formats and option combinations.

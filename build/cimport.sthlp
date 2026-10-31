{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "cimport##syntax"}{...}
{viewerjumpto "Description" "cimport##description"}{...}
{viewerjumpto "Options" "cimport##options"}{...}
{viewerjumpto "Excel options" "cimport##exceloptions"}{...}
{viewerjumpto "Other formats" "cimport##formats"}{...}
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
{synopt:{opt delimiters(chars[, collapse|asstring])}}specify field delimiter; default is automatic detection{p_end}

{syntab:Variable names}
{synopt:{opt varn:ames(rule)}}header row number or {opt nonames}; default is automatic detection{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}

{syntab:Variable types}
{synopt:{opt asfloat}}use float for noninteger numeric columns{p_end}
{synopt:{opt asdouble}}import all numeric variables as double{p_end}
{synopt:{opt numeric:cols(numlist)}}force specified columns to be numeric{p_end}
{synopt:{opt string:cols(numlist)}}force specified columns to be string; {cmd:_all} selects all columns{p_end}
{synopt:{opt favorstrfixed}}prefer fixed-width strings when values fit in 2,045 bytes{p_end}

{syntab:Parsing}
{synopt:{opt bindq:uotes(option)}}quote binding rule; {opt strict}, {opt loose}, or {opt nobind}{p_end}
{synopt:{opt stripq:uotes(default|yes|no)}}control removal of surrounding quotes{p_end}
{synopt:{opt enc:oding(encoding)}}file encoding; default is automatic detection; UTF-8/16/32 and named system encodings{p_end}
{synopt:{opt rowr:ange([start][:end])}}range of rows to import{p_end}
{synopt:{opt colr:ange([start][:end])}}range of columns to import{p_end}
{synopt:{opt empty:lines(option)}}empty line handling; {opt skip} or {opt fill}{p_end}

{syntab:Number formats}
{synopt:{opt decimalsep:arator(char)}}decimal point character; default is period{p_end}
{synopt:{opt groupsep:arator(char)}}thousands grouping character; default is none{p_end}
{synopt:{opt parsel:ocale(name)}}locale for numeric parsing{p_end}
{synopt:{opt maxquoted:rows(#)}}limit on lines in a strict quoted field{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

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
{synopt:{opt allstring}[{cmd:(}{it:format}{cmd:)}]}import all columns as string, optionally using a numeric format{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt desc:ribe}}list worksheet names and ranges without loading data{p_end}
{synopt:{opt v:erbose}}display progress information{p_end}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cimport} reads delimited text, Excel (.xls and .xlsx), dBase, SAS7BDAT,
SAS XPORT5/8, SPSS SAV/ZSAV, and ESRI shapefiles into Stata. It uses C readers
for file parsing. Text longer than 2,045 bytes is stored as {cmd:strL}.


{marker options}{...}
{title:Options for delimited text}

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

{pmore}
With several characters, each character is a delimiter. The suboption
{opt collapse} treats consecutive delimiters as one; {opt asstring} treats
the entire string as one delimiter. These suboptions cannot be combined.
For example, {cmd:delimiters("::", asstring)} splits at a pair of colons.
The top-level {opt collapsedelimiters} option is accepted but has no effect;
use {cmd:delimiters(chars, collapse)} to collapse runs.

{dlgtab:Variable names}

{phang}
{opt varnames(rule)} specifies how variable names are determined.
By default, headers are inferred. {opt varnames(1)} explicitly treats the first row as variable names.
{opt varnames(nonames)} treats the first row as data and generates default
variable names (v1, v2, ...). You can also specify {cmd:varnames(}{it:N}{cmd:)}
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
type, regardless of whether a smaller type would suffice. This retains double precision for fractional values. Cannot be combined with {opt asfloat}.

{phang}
{opt numericcols(numlist)} forces the specified columns to be imported as
numeric. Column numbers are 1-based. Values that cannot be parsed as numbers
become missing. This overrides automatic type detection for these columns.

{phang}
{opt stringcols(numlist)} forces the specified columns to be imported as
string, even if they contain only numeric values. Column numbers are 1-based.
Use {cmd:stringcols(_all)} for all columns. This overrides automatic type
detection for these columns.

{phang}
{opt favorstrfixed} uses fixed-width strings whenever the longest value fits
in 2,045 bytes. Without it, columns with highly variable text lengths may use
{cmd:strL} even when their longest value is shorter.

{dlgtab:Parsing}

{phang}
{opt bindquotes(option)} specifies how quoted fields are handled.
{opt loose} (the default) binds quoted fields within a line, but every newline
ends a row. {opt strict} allows quoted fields to span lines. {opt nobind}
ignores quotes when locating field and row boundaries.

{phang}
{opt stripquotes(default)} removes enclosing quotes when they bind a field.
{opt stripquotes(yes)} removes all quotation marks; {opt stripquotes(no)} retains them.

{phang}
{opt encoding(encoding)} overrides automatic encoding detection. Supported encodings include UTF-8,
UTF-16LE/BE, UTF-32LE/BE, ASCII, Latin-1/9, Windows code pages and Mac Roman.
{opt charset()} is an alias for {opt encoding()}.

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
Use {opt parselocale()} for delimited text; {opt locale()} is not accepted.

{phang}
{opt maxquotedrows(#)} limits physical rows within a strict quoted field and
accepts a nonnegative integer or {cmd:unlimited}. The default is 20; zero
and {cmd:unlimited} remove the limit. Exceeding the limit returns error 5101.

{dlgtab:Data recovery}

{pstd}
For delimited and Excel imports, an error or interruption restores the data
previously in memory. On Stata 16 or later, recovery usually keeps an in-memory
copy until the import finishes. Earlier versions use {help preserve}.

{phang}
{opt snapshot} uses {help preserve} to copy the data before importing. Stata
may write the copy to disk when it exceeds {cmd:c(max_preservemem)}. This
option also protects the original data for the other file formats, which
otherwise clear the dataset before reading the input.

{phang}
{opt fast} skips recovery and reduces memory use. An error or interruption can
leave an empty or incomplete import. It cannot be combined with {opt snapshot}.


{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.


{marker exceloptions}{...}
{title:Options for Excel}

{phang}
{opt sheet(name)} selects a worksheet. By default, the first worksheet is used.
Chart, dialog, and macro sheets contain no importable worksheet data.

{phang}
{opt cellrange(range)} selects a rectangular range, for example {cmd:A1:F100}.
A single cell such as {cmd:B5} sets the upper-left corner and reads through
the worksheet's remaining used range.

{phang}
{opt firstrow} uses the first imported row as variable names. Without it, names
are Excel column letters. Invalid header names are replaced with column letters;
header text is retained as variable labels. {opt case(preserve)} is the default;
{opt case(lower)} and {opt case(upper)} change variable-name case.

{phang}
{opt allstring} imports every column as text. {cmd:allstring(%12.2f)}, for
example, specifies how numeric cells are converted to strings.

{phang}
An Excel column list selects and names columns, for example
{cmd:cimport excel id=A amount=C using data.xlsx, clear}. It cannot be combined
with {opt firstrow} or {opt case()}.

{phang}
{opt describe} lists sheet names and used ranges without changing the data.
Use it without column selections or other import options. Empty worksheets
and non-worksheet tabs have no reported range.

{phang}
{opt clear}, {opt snapshot}, and {opt fast} have the meanings given above.
{opt verbose} is accepted. Excel {opt locale()} and {opt detail} are accepted
for compatibility and do not change the import. Excel import does not accept
{opt threads()}; it uses the {help ctools##threads:ctools thread setting}.

{marker formats}{...}
{title:Other file formats}

{p 8 12 2}
{cmd:cimport} {it:format} {cmd:using} {it:filename}
[{cmd:,} {opt clear} {opt snapshot} {opt fast} {it:format_options}]

{pstd}
The following formats and options are available:

{synoptset 22}{...}
{synopt:{cmd:dbase}}.dbf files; {opt case(preserve|lower|upper)}{p_end}
{synopt:{cmd:sas}}.sas7bdat files; {opt bcat(filename)}, {opt encoding(name)}, {opt case()}{p_end}
{synopt:{cmd:spss}}.sav or .zsav files; {opt zsav}, {opt encoding(name)}, {opt case()}{p_end}
{synopt:{cmd:sasxport5}}.xpt files; {opt member(name)}, {opt novallabels}, {opt describe}{p_end}
{synopt:{cmd:sasxport8}}.v8xpt files; {opt vlabfile(filename)}, {opt case()}{p_end}
{synopt:{cmd:shp}}.shp geometry files{p_end}

{pstd}
{cmd:sas} and {cmd:spss} also accept a variable list before {cmd:using} and
{cmd:if}/{cmd:in} selectors evaluated on the imported data. Other formats
in this table import the whole file. Variable-name case defaults to
{cmd:preserve}, except XPORT5, which uses lowercase names.

{pstd}
{opt bcat(filename)} supplies a SAS7BCAT catalog containing value labels.
{opt encoding(name)} supplies the source character encoding.
{opt zsav} selects the .zsav extension when the SPSS filename has no extension.

{pstd}
For XPORT5, {opt member(name)} selects a member of the transport file.
{opt describe} lists members and variables without importing; it cannot be
combined with {opt clear} or {opt novallabels}. By default, value labels are
read from {cmd:formats.xpf} in the input directory. {opt novallabels} skips
that companion file.

{pstd}
For XPORT8, {opt vlabfile(filename)} reads a companion transport file first.
If it contains data, the command stops with error 4 and leaves those data in
memory; use {opt snapshot} to retain the original dataset on this error.

{pstd}
Shapefile import reads geometry into {cmd:_ID}, {cmd:_X}, {cmd:_Y},
{cmd:rec_header}, and, for shapes that use them, {cmd:_Z} and {cmd:_M}.
The dataset retains shape metadata needed by {help cexport}. Import associated
dBase attributes separately with {cmd:cimport dbase}.

{marker remarks}{...}
{title:Excel values and dates}

{pstd}
Excel dates are converted from the workbook's 1900 or 1904 date system to
Stata dates. A numeric column containing only date-formatted temporal cells
uses daily dates; if any temporal cell includes a time, it uses datetime
values in milliseconds. All numeric cells in that column, including cells
with General formatting, are converted to the same unit. In the 1900 system,
serial 60 maps to 28 February 1900.

{pstd}
When all numeric cells share a format, the corresponding Stata display format
is retained where supported. Mixed date formats use {cmd:%td} or {cmd:%tc};
other unsupported numeric formats use {cmd:%10.0g}. With {opt allstring}, date
and time cells use their formatted Stata values, truncated to whole seconds.

{pstd}
Boolean cells become 1 or 0, except in date/datetime columns, where they become
missing. Error cells become missing. Formula cells use the results saved in
the workbook; the importer does not recalculate formulas. Unicode text and
cell text longer than 32,766 bytes are retained.


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


{pstd}
Excel and the other file formats return {cmd:r(N)}, {cmd:r(k)},
{cmd:r(width)} (bytes per observation), {cmd:r(changed)}, {cmd:r(datalabel)},
and {cmd:r(filename)}. They do not return {cmd:r(time)}.

{pstd}
{cmd:cimport excel, describe} returns {cmd:r(N_worksheet)} and the macros
{cmd:r(worksheet_1)}, {cmd:r(range_1)}, and corresponding names and ranges for
subsequent tabs. XPORT5 {opt describe} returns {cmd:r(N)}, {cmd:r(k)}, and
{cmd:r(size)} for the last member listed; without {opt member()}, it also returns
{cmd:r(n_members)} and {cmd:r(members)}.

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

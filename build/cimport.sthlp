{smcl}
{* *! version 1.2.6 20261111}{...}
{viewerjumpto "Syntax" "cimport##syntax"}{...}
{viewerjumpto "Description" "cimport##description"}{...}
{viewerjumpto "Options" "cimport##options"}{...}
{viewerjumpto "Excel options" "cimport##exceloptions"}{...}
{viewerjumpto "Parquet" "cimport##parquet"}{...}
{viewerjumpto "Other formats" "cimport##formats"}{...}
{viewerjumpto "FRED and Haver" "cimport##fred"}{...}
{viewerjumpto "Remarks" "cimport##remarks"}{...}
{viewerjumpto "Examples" "cimport##examples"}{...}
{viewerjumpto "Stored results" "cimport##results"}{...}
{title:Title}

{phang}
{bf:cimport} {hline 2} C-accelerated text, Excel, Parquet, and statistical file import


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
{synopt:{opt asdouble}}import noninteger numeric variables as double{p_end}
{synopt:{opt numeric:cols(numlist)}}force specified columns to be numeric{p_end}
{synopt:{opt string:cols(numlist)}}force specified columns to be string; {cmd:_all} selects all columns{p_end}
{synopt:{opt favorstrfixed}}prefer fixed-width strings when values fit in 2,045 bytes{p_end}

{syntab:Parsing}
{synopt:{opt bindq:uotes(option)}}quote binding rule; {opt strict}, {opt loose}, or {opt nobind}{p_end}
{synopt:{opt stripq:uotes(default|yes|no)}}control removal of surrounding quotes{p_end}
{synopt:{opt enc:oding(encoding)}}file encoding; default is automatic detection; UTF-8/16/32 and named system encodings{p_end}
{synopt:{opt rowr:ange([start][:end])}}range of rows to import{p_end}
{synopt:{opt colr:ange([start][:end])}}range of columns to import{p_end}
{synopt:{opt emptyl:ines(option)}}empty line handling; {opt skip} or {opt include}{p_end}

{syntab:Number formats}
{synopt:{opt decimals:eparator(char)}}decimal point character; default is period{p_end}
{synopt:{opt groups:eparator(char)}}thousands grouping character; default is none{p_end}
{synopt:{opt parsel:ocale(name)}}locale for numeric parsing{p_end}
{synopt:{opt maxquoted:rows(#)}}limit on lines in a strict quoted field{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

{syntab:Speed}
{synopt:{opt tempdta}}load strL-heavy data through a temporary .dta file{p_end}

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
{synopt:{opt sh:eet(name)}}worksheet name to import{p_end}
{synopt:{opt cellr:ange(range)}}cell range to import (e.g., A1:F100){p_end}
{synopt:{opt first:row}}treat first row as variable names{p_end}
{synopt:{opt allstring}[{cmd:(}{it:format}{cmd:)}]}import all columns as string, optionally using a numeric format{p_end}
{synopt:{opt case(option)}}variable name case; {opt preserve}, {opt lower}, or {opt upper}{p_end}
{synopt:{opt clear}}clear data in memory before loading{p_end}
{synopt:{opt desc:ribe}}list worksheet names and ranges without loading data{p_end}
{synopt:{opt v:erbose}}accepted; Excel import displays no timing details{p_end}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}
{synoptline}

{pstd}Import or describe a Parquet file (see {help cimport##parquet:Parquet}):

{p 8 17 2}
{cmdab:cimport}
{cmd:parquet}
[{cmd:using}]
{it:filename}
[{cmd:,} {it:parquet_options}]

{p 8 17 2}
{cmdab:cimport}
{cmd:parquet}
{it:columnlist}
{cmd:using}
{it:filename}
[{cmd:,} {it:parquet_options}]

{p 8 17 2}
{cmdab:cimport}
{cmd:parquet}
[{cmd:using}]
{it:filename}{cmd:,}
{opt desc:ribe}

{synoptset 32 tabbed}{...}
{synopthdr:parquet options}
{synoptline}
{synopt:{opt clear}}replace data in memory{p_end}
{synopt:{opt rowr:ange([start][:end])}}row range of data to load{p_end}
{synopt:{opt favormemory}}accepted for compatibility; has no effect{p_end}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}
{synoptline}

{pstd}Import from FRED or Haver Analytics with Stata's own commands (see
{help cimport##fred:FRED and Haver}):

{p 8 17 2}
{cmdab:cimport}
{cmd:fred}|{cmd:haver}|{cmd:haverdirect}
{it:...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cimport} reads delimited text, Excel (.xls, .xlsx, .xlsm, .xltx; as with
{cmd:import excel}, the kind of workbook comes from the file's content, whatever
its name), Parquet, dBase,
SAS7BDAT, SAS XPORT5/8, SPSS SAV/ZSAV, and ESRI shapefiles into Stata. It uses C
readers for file parsing. Text longer than 2,045 bytes is stored as
{cmd:strL}. {cmd:cimport fred}, {cmd:cimport haver} and {cmd:cimport haverdirect}
run Stata's {help import fred} and {help import haver}.


{marker options}{...}
{title:Options for delimited text}

{pstd}
As with {cmd:import delimited} and {cmd:import excel}, a delimited or Excel
file may be a URL ({cmd:http://}, {cmd:https://} or {cmd:ftp://}); it is copied
to a temporary file and read from there ({cmd:r(filename)} is the URL). The
other formats take local files only. A file name that starts with {cmd:~} names
a file in your home directory, as it does for Stata's own file commands.

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
{opt collapse} treats consecutive delimiters as one; {opt asstring} treats the
entire string as one delimiter. These suboptions cannot be combined. For
example, {cmd:delimiters("::", asstring)} splits at a pair of colons. The
top-level {opt collapsedelimiters} option is accepted but has no effect; use
{cmd:delimiters(chars, collapse)} to treat consecutive delimiters as one.

{dlgtab:Variable names}

{phang}
{opt varnames(rule)} specifies how variable names are determined.
By default, headers are inferred. {opt varnames(1)} explicitly treats the first row as variable names.
{opt varnames(nonames)} treats the first row as data and generates default
variable names (v1, v2, ...). You can also specify {cmd:varnames(}{it:N}{cmd:)}
to use row {it:N} as variable names (e.g., {opt varnames(3)} uses the third row).
As in {cmd:import delimited}, a name keeps the header's Unicode letters,
decimal digits and underscores, drops leading digits and is cut to 32
characters; a header that leaves no valid name gives {cmd:v}{it:#}.

{phang}
{opt case(option)} specifies the case of variable names. Delimited import defaults
to {opt lower}; Excel defaults to {opt preserve}. {opt lower} converts to lowercase.
{opt upper} converts to uppercase. Both convert Unicode letters, after the
32-character cut; a converted name longer than 32 characters falls back to
{cmd:v}{it:#} ({cmd:V}{it:#} under {opt upper}).

{dlgtab:Variable types}

{phang}
{opt asfloat} uses Stata {help data types:float} for noninteger numeric
columns. Integer columns retain their inferred byte, int, or long storage. This uses less memory than
double but has less precision. Cannot be combined with {opt asdouble}.

{phang}
{opt asdouble} uses Stata {help data types:double} for noninteger numeric
columns, as {cmd:import delimited} does. Integer columns retain their inferred
byte, int, or long storage. This retains double precision for fractional
values. Cannot be combined with {opt asfloat}.

{phang}
{opt numericcols(numlist)} forces the specified columns to be imported as
numeric. Column numbers are 1-based. Values that cannot be parsed as numbers
become missing. This overrides automatic type detection for these columns,
and, as in {cmd:import delimited}, it takes precedence over {opt stringcols()}
for a column named in both. A column holding text only is stored as
{cmd:import delimited} stores it, by its string width (1 {cmd:byte}, 2
{cmd:int}, 3 {cmd:long}, 4 {cmd:float}, 5 {cmd:double}, wider {cmd:byte}); a
column with numbers takes the type its numbers need.

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
As in {cmd:import delimited}, only a quote that is the field's first character
opens it, so {cmd: "x"} after a delimiter and a space is read as typed.
{opt stripquotes(yes)} removes all quotation marks; {opt stripquotes(no)} retains them.
Under {cmd:bindquotes(strict)}, a CRLF line ending inside quotes is kept as CR.

{pmore}
{cmd:NaN} and {cmd:Infinity}, in that case and with an optional sign, are
numeric missing values, and {cmd:d} or {cmd:D} is accepted as the exponent
marker ({cmd:1d3} is 1000), as in {cmd:import delimited}; {cmd:nan} and
{cmd:NA} are text. Numeric display formats are selected using the same sampled
lines as {cmd:import delimited}: the first 513 lines, blocks of 512 lines
beginning at lines 1,026, 2,050, 4,098, and so on, and the final incomplete
block. The format corresponds to the numeric storage type required by the
sampled values. Rows of any width are read; more columns than {cmd:maxvar}
allows fail with error 900, as in {cmd:import delimited}.

{phang}
{opt encoding(encoding)} overrides automatic encoding detection. As in
{cmd:import delimited}, it takes the names and aliases of the charsets of Java
(StataNow's runtime), in any case: {cmd:latin1}, {cmd:cyrillic},
{cmd:windows-1251}, {cmd:cp1252}, {cmd:KOI8-U}, {cmd:x-MacCyrillic}, {cmd:IBM037},
{cmd:windows-31j}, {cmd:Shift_JIS}, {cmd:GB18030}, {cmd:UTF-16} and so on.
Single-byte charsets are decoded from Java's own tables on every platform; a byte
the charset does not map is dropped, as natively. Multi-byte charsets are decoded by
the system converter. On Windows, a few have no Windows code page and give error 198:
{cmd:Big5-HKSCS}, {cmd:x-EUC-TW}, {cmd:x-ISCII91} and the IBM EBCDIC CJK pages.
{opt charset()} is an alias for {opt encoding()}.

{phang}
{opt rowrange([start][:end])} specifies a range of rows to import. As in
{cmd:import delimited}, rows are lines of the file: the variable-names line,
blank lines, and lines inside quoted fields all count. Lines end at LF, CRLF,
or a lone CR. Use
{opt rowrange(100:200)} to import lines 100-200, and {opt rowrange(100:)} or
{opt rowrange(100)} to import from line 100 to the end.

{phang}
{opt colrange([start][:end])} specifies a range of columns to import.
Use {opt colrange(2:5)} to import columns 2-5, or {opt colrange(3:)}
to import from column 3 to the last column. As in {cmd:import delimited}, the
other columns play no part: the header is detected, generated names are
numbered, and {opt numericcols()} and {opt stringcols()} count columns within
the range.

{phang}
{opt emptylines(option)} specifies how empty lines in the file are handled:
{opt skip} (the default) ignores them, and {opt include} imports each one as an
observation with all values missing, as in {cmd:import delimited}. As there,
{opt include} follows the version of the caller: under {help version} 16 or
earlier, {cmd:import delimited} reads such a file as with {opt skip}, and so
does {cmd:cimport}.

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
{opt maxquotedrows(#)} limits the number of file lines that one quoted field
may span with {cmd:bindquotes(strict)}. It accepts an integer or
{cmd:unlimited}. The default is 20; zero, a negative number and {cmd:unlimited}
remove the limit, as in {cmd:import delimited}.
Exceeding the limit returns error 5101.

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

{dlgtab:Speed}

{phang}
Without {opt tempdta}, strL values are stored through Mata and, from Stata 16,
nothing is written to disk. {opt tempdta} speeds up files whose strL text
outweighs their other data, such
as job postings with long descriptions. Stata stores strL values through Mata
several times slower than {help use} reads them, so with this option the plugin
writes the imported dataset (or, failing that, its strL columns) to a temporary
.dta file that {cmd:use} loads; the file is then erased. The file is about the
size of the imported data. It is marked temporary, which asks Windows to keep
it in memory, but it may still be written to disk, so the option is off by
default; repeated imports of large files add writes to the drive holding
{cmd:c(tmpdir)}. The result is the same as without the option, including
{cmd:c(filename)}, {cmd:c(filedate)} and {cmd:c(changed)}. The option is ignored
when the file has no strL column, when the other variables' data outweigh the
strL text (copying them through a file is slower than storing them directly),
when a column is named {cmd:_all} or {cmd:_skip}, and when an empty dataset
kept without {opt clear} holds value labels, a data label or characteristics.


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
the worksheet's remaining used range. As with {cmd:import excel}, a range that
extends beyond the worksheet's data area, a corner that is not a cell, or
corners in reverse order is an error. The used range includes cells without a
value (styled blanks, error cells, formulas without a saved result).

{phang}
{opt firstrow} uses the first imported row as variable names. Without it, names
are Excel column letters. As in {cmd:import excel}, a name keeps the Unicode
letters, decimal digits and underscores of the header's first 32 characters
and drops leading digits; a header that leaves no valid name in the requested
case is replaced with its column letter. Header text is retained as variable
labels; a numeric header cell's label is its {cmd:%18.0g} text.
{opt case(preserve)} is the default; {opt case(lower)} and {opt case(upper)}
change variable-name case, including Unicode letters.

{phang}
{opt allstring} imports every column as text ({cmd:all} is enough, as in
{cmd:import excel}). {cmd:allstring(%12.2f)}, for example, specifies how numeric
cells are converted to strings; as natively, the format is read only with the full
option name ({cmd:all(%12.2f)} imports as neither). Otherwise a number in a text
column reads as its {cmd:%18.0g} text, as in {cmd:import excel}.

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
{opt verbose} is accepted. Excel {opt locale()} is accepted and needs no
effect: Stata uses it to translate extended ASCII text (and ignores it on
Windows), while {cmd:cimport} reads text as Unicode directly, from the
workbook's own code page for .xls files. {opt detail} is accepted: the
{cmd:r(N)} and {cmd:r(k)} it adds to {cmd:import excel}'s results are always
returned. Excel import does not accept
{opt threads()}; it uses the {help ctools##threads:ctools thread setting}.

{marker parquet}{...}
{title:Parquet}

{pstd}
{cmd:cimport parquet} reads Apache Parquet files with a reader written in C,
without Java or the Arrow library. Its syntax, variables, messages and stored
results are those of StataNow's {cmd:import parquet}: the same variable names,
storage types, display formats, variable and value labels, and values.
{it:filename} without an extension means {it:filename}{cmd:.parquet}; it may
not be a URL.

{pstd}
{it:columnlist} names the columns to import, exactly as they are named in the
file; an unknown name is noted and skipped, and a name given twice is read
once. Without {opt clear}, data in memory that have changed since they were
last saved are an error. {opt describe} lists the columns and their types and
stores the numbers of rows and columns; like {cmd:import parquet}, it reads
the data and reports -1 rows for a file whose data cannot be read.

{pstd}
{opt rowrange(start:end)} loads rows {it:start} to {it:end}. {it:start} may be
{cmd:f} (the first row) and {it:end} may be {cmd:l} (the last row) or omitted;
a negative number counts from the end, so {cmd:rowrange(-10)} loads the last
10 rows.

{pstd}
Column types are converted as {cmd:import parquet} converts them:

{p2colset 9 40 42 2}{...}
{p2col:Parquet type}Stata variable{p_end}
{p2line}
{p2col:boolean}{cmd:byte} 0/1{p_end}
{p2col:int8, uint8}{cmd:byte}, or {cmd:int} when a value in the row range exceeds 100{p_end}
{p2col:int16, uint16}{cmd:int}, or {cmd:long} past 32,740{p_end}
{p2col:int32, uint32}{cmd:long}, or {cmd:double} past 2,147,483,620{p_end}
{p2col:int64, uint64}{cmd:long}; when a value in the row range lies outside
{cmd:long}'s range, {cmd:str16} holding each value's 64 bits in lowercase
hexadecimal{p_end}
{p2col:float, float16}{cmd:float}, or {cmd:double} past {cmd:c(maxfloat)}{p_end}
{p2col:double}{cmd:double}{p_end}
{p2col:date}{cmd:long}, {cmd:%td}{p_end}
{p2col:time}{cmd:double}, milliseconds since midnight (finer units truncated){p_end}
{p2col:timestamp, INT96}{cmd:double}, {cmd:%tc}, milliseconds (finer units
truncated; UTC adjustment is ignored){p_end}
{p2col:decimal}{cmd:str}{it:#} text, {it:#} the precision plus 2{p_end}
{p2col:string, JSON}{cmd:str}{it:#}, {it:#} the longest value in the file
(in bytes), or {cmd:strL} past 2,045{p_end}
{p2col:binary, enum, BSON, fixed-length}{cmd:strL}{p_end}
{p2col:Arrow dictionary}{cmd:byte}, {cmd:int} or {cmd:long} codes 0, 1, ...,
labeled with the dictionary's entries{p_end}
{p2line}

{pstd}
A value a type cannot hold is missing. Null values are missing ({cmd:.} or
empty text). Nested columns (structs, lists, maps, repeated fields), null-type
columns and Arrow durations are skipped with a note. Column names become valid
Stata names: characters other than letters, digits and underscores are
removed, leading digits are dropped and names are cut to 32 characters; an
empty name becomes {cmd:v}{it:#}, {it:#} the column's position counting from
0, and if two names would be the same, every variable is named that way. A
variable whose name differs from its column is labeled with the column's
name.

{pstd}
All row groups, data page versions 1 and 2, and the PLAIN, dictionary, RLE,
BIT_PACKED, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY and
BYTE_STREAM_SPLIT encodings are read, uncompressed or compressed with Snappy,
Gzip, Brotli, Zstandard, LZ4 or LZ4_RAW. Neither command reads LZO.
{opt favormemory} changes nothing: {cmd:cimport parquet} always decodes only
the selected columns.

{pstd}
By default, an import that fails restores the data in memory; {opt fast}
skips that copy, and {opt snapshot} makes it with {help preserve}.


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
{cmd:sas} and {cmd:spss} also accept a list of variable names and
{cmd:if}/{cmd:in} before {cmd:using}. As in {cmd:import sas} and {cmd:import spss},
the names are the file's exact (case-sensitive) names, a hyphen separates names
as a blank does ({cmd:a-b} is {cmd:a} and {cmd:b}), names the file lacks are
ignored, and with no match the import fails ("no matching names found in data",
r(111)). The names and an {cmd:in} range are applied while the file is read, so
only the selection is held in memory; {cmd:if} is then applied to the selected
variables. Other formats
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

{marker fred}{...}
{title:FRED and Haver}

{pstd}
{cmd:cimport fred}, {cmd:cimport haver} and {cmd:cimport haverdirect} (and
their abbreviations) pass all their arguments, unchanged, to Stata's
{help import fred}, {help import haver} and {help import haverdirect}. FRED is
a web service and Haver databases need Haver Analytics' Windows library, so
ctools has no C reader for them; the results, messages, stored results and
the data left in memory after an error are those of Stata's commands.
Options {opt snapshot} and {opt fast} do not apply.


{marker remarks}{...}
{title:Excel values and dates}

{pstd}
Excel dates are converted from the workbook's 1900 or 1904 date system to
Stata dates. If all date or time formats in a numeric column specify dates
only, the column uses Stata daily dates. If any specifies a time, the column
uses Stata datetime values in milliseconds. All numeric cells in that column,
including cells with General formatting, are converted to the same unit. In
the 1900 system, serial 60 maps to 28 February 1900.

{pstd}
When all numeric cells share a Stata display format, it is retained where
supported: built-in date and time formats map to their Stata formats, and a
custom format maps to {cmd:%td} (a date), {cmd:%tc} (a date with a time, or an
elapsed time such as {cmd:[h]:mm}) or {cmd:%tcHH:MM:SS} (a time of day). Mixed
date formats use {cmd:%td} or {cmd:%tc}, or {cmd:%tcHH:MM:SS} when every date
or time cell is a time; other unsupported numeric formats use {cmd:%10.0g}. A
cell without a style takes the workbook's first cell style. With
{opt allstring}, and in text columns and header labels, date and time cells use
their formatted Stata values, truncated to whole seconds; a value Stata cannot
show as a date (beyond year 9999) reads as Stata displays it.

{pstd}
Boolean cells become 1 or 0, except in date/datetime columns, where they become
missing. Error cells become missing. Formula cells use the results saved in
the workbook; the importer does not recalculate formulas. Unicode text and
cell text longer than 32,766 bytes are retained. A cell value that is not a
number reads as 0, and one that begins with a number as that number, as in
{cmd:import excel}. Values cached inside a chart embedded in an .xls worksheet
are not read as cells. In .xlsx text, a carriage return, or a carriage return
and line feed, written as such reads as a line feed, as XML requires and as
{cmd:import excel} reads it; one written as {cmd:_x000D_} or {cmd:&#13;} stays
a carriage return.

{pstd}
An .xlsx file contains compressed files. Before extracting each one,
{cmd:cimport} checks its reported size. An impossible reported size causes
error 603, as in {cmd:import excel}. A file larger than half of physical
memory, or than {cmd:c(max_memory)} when that is set, causes error 909. In
either case, the data in memory are kept.


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

{pstd}Import column 1 as a string:{p_end}
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

{pstd}{ul:Parquet import examples}

{pstd}Describe the columns of a Parquet file:{p_end}
{phang2}{cmd:. cimport parquet using auto, describe}{p_end}

{pstd}Import rows 64 through 74 of columns make and price:{p_end}
{phang2}{cmd:. cimport parquet make price using auto, rowrange(64:74) clear}{p_end}


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
Excel and the other file formats return {cmd:r(N)}, {cmd:r(k)}, {cmd:r(width)}
(bytes per observation), {cmd:r(changed)}, {cmd:r(datalabel)}, and
{cmd:r(filename)}. They do not return {cmd:r(time)}.

{pstd}
{cmd:cimport parquet} returns {cmd:r(N)} and {cmd:r(k)}, and
{cmd:cimport parquet, describe} returns {cmd:r(N_rows)} and
{cmd:r(n_columns)}, as {cmd:import parquet} does. {cmd:cimport fred} and
{cmd:cimport haver} return what Stata's commands return.

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
Manual: {bf:[D] import delimited}, {bf:[D] import excel}, {bf:[D] import parquet}

{psee}
Online: {help import delimited}, {help import excel}, {help import parquet}, {help insheet}, {help cexport}, {help ctools}
{p_end}

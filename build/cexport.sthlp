{smcl}
{* *! version 1.2.1 20261029}{...}
{viewerjumpto "Syntax" "cexport##syntax"}{...}
{viewerjumpto "Description" "cexport##description"}{...}
{viewerjumpto "Options" "cexport##options"}{...}
{viewerjumpto "Options for excel" "cexport##exceloptions"}{...}
{viewerjumpto "Other formats" "cexport##formats"}{...}
{viewerjumpto "Remarks" "cexport##remarks"}{...}
{viewerjumpto "Examples" "cexport##examples"}{...}
{viewerjumpto "Stored results" "cexport##results"}{...}
{title:Title}

{phang}
{bf:cexport} {hline 2} C-accelerated text, Excel, and statistical file export


{marker syntax}{...}
{title:Syntax}

{pstd}Export to delimited text file (CSV, TSV, etc.):

{p 8 17 2}
{cmdab:cexport}
{cmd:delimited}
[{varlist}]
{cmd:using}
{it:filename}
{ifin}
[{cmd:,} {it:options}]

{pstd}Export to Excel (.xls or .xlsx) file:

{p 8 17 2}
{cmdab:cexport}
{cmd:excel}
[{varlist}]
{cmd:using}
{it:filename}
{ifin}
[{cmd:,} {it:excel_options}]

{synoptset 24 tabbed}{...}
{synopthdr:delimited options}
{synoptline}
{syntab:Main}
{synopt:{opt d:elimiter(char)}}field delimiter; default is comma{p_end}
{synopt:{opt replace}}overwrite existing file{p_end}

{syntab:Formatting}
{synopt:{opt novarnames}}do not write variable names as header row{p_end}
{synopt:{opt quote}}quote all string fields{p_end}
{synopt:{opt noquoteif}}do not automatically quote strings containing delimiters{p_end}
{synopt:{opt nolabel}}export values instead of value labels{p_end}
{synopt:{opt datafmt}}export numeric variables using their display formats{p_end}
{synopt:{opt datestring(fmt)}}custom format for date/time variables{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress and elapsed times{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:I/O options}
{synopt:{opt mmap}}use memory-mapped I/O (zero-copy formatting){p_end}
{synopt:{opt nofsync}}skip the final disk synchronization{p_end}
{synopt:{opt direct}}request direct I/O where supported{p_end}
{synopt:{opt prefault}}pre-fault mmap pages to avoid stalls{p_end}
{synopt:{opt crlf}}use Windows-style CRLF line endings{p_end}
{synopt:{opt noparallel}}disable parallel I/O{p_end}
{synoptline}

{synoptset 24 tabbed}{...}
{synopthdr:excel options}
{synoptline}
{syntab:Main}
{synopt:{opt sheet(name[, modify|replace])}}worksheet name and update mode; default is "Sheet1"{p_end}
{synopt:{opt cell(start)}}starting cell for export; default is "A1"{p_end}
{synopt:{opt replace}}overwrite existing file{p_end}

{syntab:Formatting}
{synopt:{opt firstrow(variables|varlabels)}}write column headers; default writes data in the first row{p_end}
{synopt:{opt nolabel}}export values instead of value labels{p_end}
{synopt:{opt datafmt}}export date/time variables using their display formats{p_end}
{synopt:{opt datestring(fmt)}}custom format for date/time variables{p_end}
{synopt:{opt missing(string)}}replacement value for missing data; default is empty cell{p_end}
{synopt:{opt keepcellfmt}}retain existing cell formats when modifying a worksheet{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synopt:{opt threads(#)}}maximum number of threads to use{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cexport} writes delimited text, Excel (.xls and .xlsx), dBase, SAS7BDAT,
SAS XPORT5/8, SPSS SAV, and ESRI shapefiles. It uses C routines for file
formatting and writing. Unless a variable list or {cmd:if}/{cmd:in} is specified,
all variables and observations are exported. Shapefile export uses its own
geometry schema, described below.


{marker options}{...}
{title:Options for delimited}

{dlgtab:Main}

{phang}
{opt delimiter(char)} specifies the delimiter to use between fields. The
default is comma ({cmd:,}). Use {cmd:delimiter(tab)} or {cmd:delimiter(\t)}
for tab-delimited output.

{phang}
{opt replace} specifies that {it:filename} be replaced if it already exists.

{dlgtab:Formatting}

{phang}
{opt novarnames} specifies that variable names should not be written as the
first row of the file.

{phang}
{opt quote} specifies that all string fields should be enclosed in double
quotes, regardless of content.

{phang}
{opt noquoteif} specifies that strings should never be quoted, even if they
contain the delimiter character.

{phang}
{opt nolabel} specifies that numeric values should be exported as raw values
instead of their value labels.

{phang}
{opt datafmt} writes numeric and date/time values using their Stata display
formats. Without it, delimited export writes numeric values without applying
display formats.


{phang}
{opt datestring(fmt)} specifies a custom format to use for all date/time
variables. This overrides the variables' display formats. Common formats include:

{p 12 16 2}{cmd:datestring("%tdCCYY-NN-DD")} for ISO 8601 dates (2024-01-15){p_end}
{p 12 16 2}{cmd:datestring("%tdNN/DD/CCYY")} for US format (01/15/2024){p_end}
{p 12 16 2}{cmd:datestring("%tdDD-Mon-CCYY")} for European format (15-Jan-2024){p_end}
{p 12 16 2}{cmd:datestring("%tcCCYY-NN-DD!THH:MM:SS")} for ISO 8601 datetime{p_end}

{dlgtab:Reporting}

{phang}
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{dlgtab:I/O options}

{phang}
{opt mmap} uses memory-mapped I/O for zero-copy formatting, which can improve
performance for large files.

{phang}
{opt nofsync} skips the final fsync system call for faster writes. The file
may not be fully flushed to disk when the command returns. The operating system may still have buffered data if the machine loses power.

{phang}
{opt direct} requests direct I/O to bypass the operating system page cache. Support depends on the operating system and selected I/O backend.

{phang}
{opt prefault} pre-faults mmap pages to avoid page fault stalls during writing.

{phang}
{opt crlf} uses Windows-style CRLF ({cmd:\r\n}) line endings instead of the
default Unix-style LF ({cmd:\n}).

{phang}
{opt noparallel} disables parallel I/O.


{marker exceloptions}{...}
{title:Options for excel}

{dlgtab:Main}

{phang}
{opt sheet(name)} specifies the worksheet name. The default is "Sheet1".
The name is limited to 31 characters. With an existing workbook,
{cmd:sheet(name, modify)} writes into that worksheet while retaining cells
outside the exported rectangle. {cmd:sheet(name, replace)} replaces the entire
worksheet and retains other sheets. The top-level options {opt sheetmodify}
and {opt sheetreplace} select the same modes. Neither mode may be combined
with {opt replace}, which replaces the entire workbook.

{phang}
{opt cell(start)} specifies the starting cell for the data export. The default
is "A1". For example, {cmd:cell(B5)} starts the export at column B, row 5.
This allows you to write data to a specific region of a worksheet.

{phang}
{opt replace} specifies that {it:filename} be replaced if it already exists.

{dlgtab:Formatting}

{phang}
By default, Excel export writes data starting in the first row.
{cmd:firstrow(variables)} writes variable names as column headers.
{cmd:firstrow(varlabels)} writes variable labels, falling back to names.
{cmd:firstrow(nonames)} remains an extension for explicitly omitting headers.

{phang}
{opt nolabel} specifies that numeric values should be exported as raw values
instead of their value labels.

{phang}
{opt missing(string)} specifies a string value to use for missing data.
By default, missing values are exported as empty cells. For example,
{cmd:missing("NA")} exports all missing values as "NA", and {cmd:missing(".")}
exports them as periods.

{phang}
{opt keepcellfmt} retains existing formats for cells overwritten by
{cmd:sheet(name, modify)}. Without it, overwritten cells receive the export's
formats. It has no effect on a new workbook. With whole-workbook {opt replace},
it retains style definitions but not the old cells, their format assignments,
or other worksheets; use sheet modification to retain an existing layout.

{phang}
{opt datafmt} writes date/time variables as strings in their Stata display
formats. {opt datestring(fmt)} supplies a common date/time format instead.
These options do not format ordinary numeric variables for Excel.
Without them, dates are converted to Excel date values and assigned date
formats; see {help cexport##remarks:Remarks}.

{phang}
{opt locale(name)} is accepted for compatibility and has no effect.


{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information.


{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{marker formats}{...}
{title:Other file formats}

{p 8 12 2}
{cmd:cexport} {it:format} [{varlist}] {ifin} {cmd:using} {it:filename}
[{cmd:,} {opt replace} {it:format_options}]

{synoptset 22}{...}
{synopt:{cmd:dbase}}.dbf files; {opt version(3)}, {opt datafmt}, {opt origdbfdate}{p_end}
{synopt:{cmd:sas}}.sas7bdat files{p_end}
{synopt:{cmd:spss}}.sav files; {opt novallabels}{p_end}
{synopt:{cmd:sasxport5}}.xpt files; {opt rename}, {opt vallabfile(xpf|sascode|both|none)}{p_end}
{synopt:{cmd:sasxport8}}.v8xpt files; {opt vallabfile}{p_end}

{pstd}
{opt novallabels} omits SPSS value labels. XPORT5 uses {cmd:vallabfile(xpf)} by
default, writing labels to {cmd:formats.xpf} in the destination directory.
{cmd:sascode} writes a .sas program, {cmd:both} writes both files, and
{cmd:none} writes neither. {opt rename} converts names to XPORT5's eight-character
limit. XPORT8 {opt vallabfile} writes a .sas label program.

{pstd}
dBase export defaults to version IV; {cmd:version(3)} selects version III.
{opt datafmt} applies numeric display formats. {opt origdbfdate} uses the date
retained in dataset characteristics by a prior dBase import.

{p 8 12 2}
{cmd:cexport shp using} {it:filename}
[{cmd:,} {opt replace} {opt shx} {opt id(varname)}]

{pstd}
Shapefile export requires the geometry variables and dataset characteristics
from a shapefile import: {cmd:_ID}, {cmd:_X}, {cmd:_Y}, {cmd:rec_header}, and
{cmd:_Z}/{cmd:_M} where used by the shape type. {opt id()} names an alternative
numeric identifier. {opt shx} also writes the .shx index. This syntax exports
the whole geometry dataset and does not accept a variable list or {cmd:if}/{cmd:in}.

{marker remarks}{...}
{title:Remarks}

{pstd}
Delimited and Excel exports write a temporary file and publish the destination
after successful completion. Without {opt replace} or a sheet-update option,
an existing destination is an error.

{pstd}
New Excel workbooks use the 1900 date system. Daily dates retain their calendar
date; datetimes retain fractional days, including milliseconds. {cmd:%tC}
values are converted with {cmd:cofC()} because Excel has no leap-second
representation. Other Stata period dates use the first day of the period.

{pstd}
dBase numeric fields use widths of 3, 5, and 10 for {cmd:byte}, {cmd:int}, and
{cmd:long}, widened by one when a negative value needs it. Other numeric
variables use 20 columns with {cmd:%20.0g} formatting unless {opt datafmt}
is specified.

{marker examples}{...}
{title:Examples}

{pstd}Setup:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}{ul:Delimited export examples}

{pstd}Export all variables to a CSV file:{p_end}
{phang2}{cmd:. cexport delimited using auto.csv, replace}{p_end}

{pstd}Export selected variables:{p_end}
{phang2}{cmd:. cexport delimited make price mpg using auto_subset.csv, replace}{p_end}

{pstd}Export to tab-delimited file:{p_end}
{phang2}{cmd:. cexport delimited using auto.tsv, delimiter(tab) replace}{p_end}

{pstd}Export without header row:{p_end}
{phang2}{cmd:. cexport delimited using auto_noheader.csv, novarnames replace}{p_end}

{pstd}Export with verbose output:{p_end}
{phang2}{cmd:. cexport delimited using auto.csv, replace verbose}{p_end}

{pstd}Export subset of observations:{p_end}
{phang2}{cmd:. cexport delimited using foreign_cars.csv if foreign == 1, replace}{p_end}

{pstd}Export with date formatting (using variable display formats):{p_end}
{phang2}{cmd:. cexport delimited using auto.csv, datafmt replace}{p_end}

{pstd}{ul:Excel export examples}

{pstd}Export all variables to an Excel file:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, replace}{p_end}

{pstd}Export with custom sheet name:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, sheet("AutoData") replace}{p_end}

{pstd}Export variable names in the first row:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, firstrow(variables) replace}{p_end}

{pstd}Update cells in an existing worksheet:{p_end}
{phang2}{cmd:. cexport excel price mpg using auto.xlsx, sheet("Sheet1", modify) cell(B2) keepcellfmt}{p_end}

{pstd}Export selected variables:{p_end}
{phang2}{cmd:. cexport excel make price mpg weight using auto_subset.xlsx, replace}{p_end}

{pstd}Export subset of observations:{p_end}
{phang2}{cmd:. cexport excel using domestic.xlsx if foreign == 0, replace}{p_end}

{pstd}Export starting at cell B5:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, cell(B5) replace}{p_end}

{pstd}Export with custom missing value:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, missing("NA") replace}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cexport delimited} and {cmd:cexport excel} store the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations exported{p_end}
{synopt:{cmd:r(k)}}number of variables exported{p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(filename)}}name of the output file{p_end}


{pstd}
Other statistical and dBase exports return {cmd:r(N)}, {cmd:r(k)}, and
{cmd:r(filename)}. Shapefile export returns {cmd:r(filename)} only.

{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] export delimited}, {bf:[D] export excel}

{psee}
Online: {help export delimited}, {help export excel}, {help cimport}, {help ctools}
{p_end}

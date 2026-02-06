{smcl}
{* *! version 0.9.1 26Jan2026}{...}
{viewerjumpto "Syntax" "cexport##syntax"}{...}
{viewerjumpto "Description" "cexport##description"}{...}
{viewerjumpto "Options" "cexport##options"}{...}
{viewerjumpto "Options for excel" "cexport##exceloptions"}{...}
{viewerjumpto "Examples" "cexport##examples"}{...}
{viewerjumpto "Stored results" "cexport##results"}{...}
{title:Title}

{phang}
{bf:cexport} {hline 2} C-accelerated data export (CSV and Excel)


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

{pstd}Export to Excel (.xlsx) file:

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
{synopt:{opt datafmt}}export date/time variables using their display formats{p_end}
{synopt:{opt datestring(fmt)}}custom format for date/time variables{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synopt:{opt timeit}}display timing breakdown{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}

{syntab:Advanced Performance}
{synopt:{opt mmap}}use memory-mapped I/O (zero-copy formatting){p_end}
{synopt:{opt nofsync}}skip final fsync for faster writes (less durable){p_end}
{synopt:{opt direct}}use direct I/O bypassing OS cache (for very large files){p_end}
{synopt:{opt prefault}}pre-fault mmap pages to avoid stalls{p_end}
{synopt:{opt crlf}}use Windows-style CRLF line endings{p_end}
{synopt:{opt noparallel}}disable parallel I/O (for debugging){p_end}
{synoptline}

{synoptset 24 tabbed}{...}
{synopthdr:excel options}
{synoptline}
{syntab:Main}
{synopt:{opt sheet(name)}}worksheet name; default is "Sheet1"{p_end}
{synopt:{opt cell(start)}}starting cell for export; default is "A1"{p_end}
{synopt:{opt replace}}overwrite existing file{p_end}

{syntab:Formatting}
{synopt:{opt firstrow(variables|nonames)}}first row contains variable names (default) or data{p_end}
{synopt:{opt nolabel}}export values instead of value labels{p_end}
{synopt:{opt datafmt}}export date/time variables using their display formats{p_end}
{synopt:{opt datestring(fmt)}}custom format for date/time variables{p_end}
{synopt:{opt missing(string)}}replacement value for missing data; default is empty cell{p_end}
{synopt:{opt keepcellfmt}}preserve cell formatting from existing file{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cexport} provides high-performance data export using C plugins with parallel
processing. It supports two formats:

{phang2}
{cmd:cexport delimited} exports data to delimited text files (CSV, TSV, etc.).
It is a high-performance replacement for {help export delimited:export delimited}.

{phang2}
{cmd:cexport excel} exports data to Excel (.xlsx) files. It is a high-performance
replacement for {help export excel:export excel}.


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
{opt datafmt} specifies that date/time variables should be exported using their
display formats (e.g., {cmd:%td}, {cmd:%tc}, {cmd:%tw}). Without this option,
date/time variables are exported as raw numeric values (days since 1960-01-01
for daily dates, milliseconds since 1960-01-01 for datetimes, etc.).

{phang}
{opt datestring(fmt)} specifies a custom format to use for all date/time
variables. This overrides the variables' display formats. Common formats include:

{p 12 16 2}{cmd:datestring("%tdCCYY-NN-DD")} for ISO 8601 dates (2024-01-15){p_end}
{p 12 16 2}{cmd:datestring("%tdNN/DD/CCYY")} for US format (01/15/2024){p_end}
{p 12 16 2}{cmd:datestring("%tdDD-Mon-CCYY")} for European format (15-Jan-2024){p_end}
{p 12 16 2}{cmd:datestring("%tcCCYY-NN-DD!THH:MM:SS")} for ISO 8601 datetime{p_end}

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information during export.

{phang}
{opt timeit} displays a timing breakdown showing data loading time, write
time, and throughput.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
data loading and export operations. By default, {cmd:cexport} uses all available
CPU cores as reported by OpenMP. Use this option to limit parallelism, for
example when running multiple jobs simultaneously.

{dlgtab:Advanced Performance}

{phang}
{opt mmap} uses memory-mapped I/O for zero-copy formatting, which can improve
performance for large files.

{phang}
{opt nofsync} skips the final fsync system call for faster writes. The file
may not be fully flushed to disk when the command returns. Use when speed
matters more than durability.

{phang}
{opt direct} uses direct I/O, bypassing the OS page cache. Useful for
very large files to avoid evicting other data from cache.

{phang}
{opt prefault} pre-faults mmap pages to avoid page fault stalls during writing.

{phang}
{opt crlf} uses Windows-style CRLF ({cmd:\r\n}) line endings instead of the
default Unix-style LF ({cmd:\n}).

{phang}
{opt noparallel} disables parallel I/O. Useful for debugging.


{marker exceloptions}{...}
{title:Options for excel}

{dlgtab:Main}

{phang}
{opt sheet(name)} specifies the worksheet name. The default is "Sheet1".
The name is limited to 31 characters.

{phang}
{opt cell(start)} specifies the starting cell for the data export. The default
is "A1". For example, {cmd:cell(B5)} starts the export at column B, row 5.
This allows you to write data to a specific region of a worksheet.

{phang}
{opt replace} specifies that {it:filename} be replaced if it already exists.

{dlgtab:Formatting}

{phang}
{opt firstrow(variables|nonames)} specifies how to handle the first row.
{cmd:firstrow(variables)} (the default) writes variable names as the first row.
{cmd:firstrow(nonames)} writes data starting from the first row.

{phang}
{opt nolabel} specifies that numeric values should be exported as raw values
instead of their value labels.

{phang}
{opt missing(string)} specifies a string value to use for missing data.
By default, missing values are exported as empty cells. For example,
{cmd:missing("NA")} exports all missing values as "NA", and {cmd:missing(".")}
exports them as periods.

{phang}
{opt keepcellfmt} preserves cell formatting (fonts, colors, borders, number
formats) from an existing Excel file when replacing it. This allows you to
update data in a formatted template without losing the formatting. The option
has no effect when creating a new file. Requires the {cmd:replace} option.

{dlgtab:Reporting}

{phang}
{opt verbose} displays timing information during export.


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

{pstd}Export with verbose timing output:{p_end}
{phang2}{cmd:. cexport delimited using auto.csv, replace verbose timeit}{p_end}

{pstd}Export subset of observations:{p_end}
{phang2}{cmd:. cexport delimited using foreign_cars.csv if foreign == 1, replace}{p_end}

{pstd}Export with date formatting (using variable display formats):{p_end}
{phang2}{cmd:. cexport delimited using auto.csv, datafmt replace}{p_end}

{pstd}{ul:Excel export examples}

{pstd}Export all variables to an Excel file:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, replace}{p_end}

{pstd}Export with custom sheet name:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, sheet("AutoData") replace}{p_end}

{pstd}Export without variable names in first row:{p_end}
{phang2}{cmd:. cexport excel using auto.xlsx, firstrow(nonames) replace}{p_end}

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
{cmd:cexport} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations exported{p_end}
{synopt:{cmd:r(k)}}number of variables exported{p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(filename)}}name of the output file{p_end}


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

{smcl}
{* *! version 1.1.0}{...}
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

{syntab:Reporting}
{synopt:{opt verbose}}display progress information{p_end}
{synopt:{opt timeit}}display timing breakdown{p_end}
{synoptline}

{synoptset 24 tabbed}{...}
{synopthdr:excel options}
{synoptline}
{syntab:Main}
{synopt:{opt sheet(name)}}worksheet name; default is "Sheet1"{p_end}
{synopt:{opt replace}}overwrite existing file{p_end}

{syntab:Formatting}
{synopt:{opt firstrow(variables|nonames)}}first row contains variable names (default) or data{p_end}
{synopt:{opt nolabel}}export values instead of value labels{p_end}

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

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information during export.

{phang}
{opt timeit} displays a timing breakdown showing data loading time, write
time, and throughput.


{marker exceloptions}{...}
{title:Options for excel}

{dlgtab:Main}

{phang}
{opt sheet(name)} specifies the worksheet name. The default is "Sheet1".
The name is limited to 31 characters.

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

{dlgtab:Reporting}

{phang}
{opt verbose} displays timing information during export.


{marker examples}{...}
{title:Examples}

{pstd}{ul:Delimited export examples}

{pstd}Export all variables to a CSV file:{p_end}
{phang2}{cmd:. cexport delimited using output.csv, replace}{p_end}

{pstd}Export selected variables:{p_end}
{phang2}{cmd:. cexport delimited id name value using output.csv, replace}{p_end}

{pstd}Export to tab-delimited file:{p_end}
{phang2}{cmd:. cexport delimited using output.tsv, delimiter(tab) replace}{p_end}

{pstd}Export without header row:{p_end}
{phang2}{cmd:. cexport delimited using output.csv, novarnames replace}{p_end}

{pstd}Export with verbose timing output:{p_end}
{phang2}{cmd:. cexport delimited using output.csv, replace verbose timeit}{p_end}

{pstd}Export subset of observations:{p_end}
{phang2}{cmd:. cexport delimited using subset.csv if year > 2020, replace}{p_end}

{pstd}{ul:Excel export examples}

{pstd}Export all variables to an Excel file:{p_end}
{phang2}{cmd:. cexport excel using output.xlsx, replace}{p_end}

{pstd}Export with custom sheet name:{p_end}
{phang2}{cmd:. cexport excel using output.xlsx, sheet("Data") replace}{p_end}

{pstd}Export without variable names in first row:{p_end}
{phang2}{cmd:. cexport excel using output.xlsx, firstrow(nonames) replace}{p_end}

{pstd}Export selected variables:{p_end}
{phang2}{cmd:. cexport excel id name value using output.xlsx, replace}{p_end}

{pstd}Export subset of observations:{p_end}
{phang2}{cmd:. cexport excel using subset.xlsx if year > 2020, replace}{p_end}


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
ctools package


{title:Also see}

{psee}
{space 2}Help: {help export delimited}, {help export excel}, {help ctools}, {help cimport}
{p_end}

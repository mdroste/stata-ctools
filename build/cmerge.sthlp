{smcl}
{* *! version 1.0.1 07Feb2026}{...}
{viewerjumpto "Syntax" "cmerge##syntax"}{...}
{viewerjumpto "Description" "cmerge##description"}{...}
{viewerjumpto "Options" "cmerge##options"}{...}
{viewerjumpto "Remarks" "cmerge##remarks"}{...}
{viewerjumpto "Examples" "cmerge##examples"}{...}
{viewerjumpto "Stored results" "cmerge##results"}{...}
{title:Title}

{phang}
{bf:cmerge} {hline 2} C-accelerated merge for Stata datasets


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cmerge}
{it:merge_type}
{varlist}
{cmd:using}
{it:filename}
[{cmd:,} {it:options}]

{pstd}
where {it:merge_type} is one of

{p2colset 9 22 24 2}{...}
{p2col :{opt 1:1}}one-to-one merge{p_end}
{p2col :{opt m:1}}many-to-one merge{p_end}
{p2col :{opt 1:m}}one-to-many merge{p_end}
{p2col :{opt m:m}}many-to-many merge (not recommended){p_end}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt keep(results)}}which observations to keep from merged data{p_end}
{synopt:{opt assert(results)}}verify that merge results are as expected{p_end}
{synopt:{opt gen:erate(varname)}}name of variable to mark merge results; default is {cmd:_merge}{p_end}
{synopt:{opt nogen:erate}}do not create {cmd:_merge} variable{p_end}
{synopt:{opt keepus:ing(varlist)}}variables to keep from using data{p_end}
{syntab:Merge behavior}
{synopt:{opt sorted}}assert that both datasets are sorted on key variables{p_end}
{synopt:{opt force}}allow string/numeric mismatches for non-key variables only{p_end}
{synopt:{opt nolabel}}do not copy value labels from using data{p_end}
{synopt:{opt nonotes}}do not copy variable notes from using data{p_end}
{synopt:{opt update}}update missing values of same-named variables with using data{p_end}
{synopt:{opt replace}}replace all values of same-named variables with using data{p_end}
{synopt:{opt preserve_order(#)}}preserve observation order from master dataset (0=no, 1=yes){p_end}
{syntab:Reporting}
{synopt:{opt norep:ort}}do not display result summary{p_end}
{synopt:{opt verbose}}display detailed timing and progress information{p_end}
{synopt:{opt time:it}}display timing breakdown{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{pstd}
{it:results} is one or more of: {opt match}, {opt master}, {opt using}, {opt 1}, {opt 2}, {opt 3}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cmerge} is a high-performance drop-in replacement for Stata's {help merge:merge}
command. It uses a C plugin with parallel data loading and optimized sorting
to achieve significant speed improvements over the native command.

{pstd}
{cmd:cmerge} joins the dataset currently in memory (the master data) with
{it:filename} (the using data), matching on the key variables specified in
{it:varlist}. The result replaces the data in memory.


{marker options}{...}
{title:Options}

{phang}
{opt keep(results)} specifies which observations are to be kept from the
merged data. Valid results are {opt match} (or {opt 3}), {opt master} (or {opt 1}),
and {opt using} (or {opt 2}).

{phang}
{opt assert(results)} specifies assertions about the match results that should
be verified. If any assertion fails, an error is raised.

{phang}
{opt generate(varname)} specifies the name of the variable to be created
marking the source of each observation. The default name is {cmd:_merge}.

{phang}
{opt nogenerate} specifies that the merge-result variable should not be created.

{phang}
{opt keepusing(varlist)} specifies which variables from the using dataset
should be kept in the merged result.

{phang}
{opt sorted} specifies that both datasets are already sorted on the key
variables. This skips the internal sorting step and can improve performance
for pre-sorted data.

{phang}
{opt force} allows incompatible string/numeric types in shared non-key variables.
The master storage type is retained, and incompatible using-side values are
missing. Key variables must have matching string/numeric types even with
{opt force}; otherwise the command returns error 106 before modifying data.
{p_end}
{p 8 12 2}2. High-performance radix sort on key variables{p_end}
{p 8 12 2}3. Single-pass sorted merge join{p_end}
{p 8 12 2}4. Parallel output construction{p_end}

{pstd}
Runtime depends on dataset shape, options, and hardware. Benchmark comparisons
should record the date, CPU/OS, Stata/reference versions, options, and dimensions.

{pstd}
Due to the overhead of reading data from Stata into C memory and writing results
back, the largest speedups occur when both the master and using datasets (or the
variables specified in {opt keepusing()} if used) have relatively few variables
and the dataset has many observations. For merges involving many variables,
the data transfer overhead may reduce the relative performance advantage.


{marker examples}{...}
{title:Examples}

{pstd}Create example datasets for merging:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. keep make price mpg}{p_end}
{phang2}{cmd:. save auto_master, replace}{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. keep make weight length foreign}{p_end}
{phang2}{cmd:. save auto_using, replace}{p_end}

{pstd}Simple one-to-one merge on make:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using}{p_end}

{pstd}Keep only matched observations:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, keep(match)}{p_end}

{pstd}Merge with verbose output and no merge variable:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, nogenerate verbose}{p_end}

{pstd}Keep only selected variables from using dataset:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, keepusing(weight foreign)}{p_end}

{pstd}Many-to-one merge example with panel data:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. keep idcode year ln_wage}{p_end}
{phang2}{cmd:. save nlswork_wages, replace}{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. collapse (first) race, by(idcode)}{p_end}
{phang2}{cmd:. save nlswork_person, replace}{p_end}
{phang2}{cmd:. use nlswork_wages, clear}{p_end}
{phang2}{cmd:. cmerge m:1 idcode using nlswork_person}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cmerge} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations after merge{p_end}
{synopt:{cmd:r(N_1)}}number of observations only in master{p_end}
{synopt:{cmd:r(N_2)}}number of observations only in using{p_end}
{synopt:{cmd:r(N_3)}}number of matched observations{p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(using)}}name of using file{p_end}
{synopt:{cmd:r(keyvars)}}key variable names{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] merge}

{psee}
Online: {help merge}, {help joinby}, {help ctools}, {help csort}
{p_end}

{title:Validation and failure behavior}

{pstd}
Shared numeric and fixed-string variables are widened before updates or using-only rows are written. A long/float combination uses double to preserve both inputs. Key variables must have matching string/numeric types even with force. For incompatible non-key types, force retains the master type and treats the using-side values as missing. String writes retain the shared plugin limit of str2045; strL output is unsupported. Failed merges restore the master dataset.
{p_end}
